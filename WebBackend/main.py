from flask import Flask, request, jsonify, send_file
from flask_cors import CORS
from TechnicalToolsV2 import log
from file_handler import *
from time import time, sleep

import SunburstChartSetup
import secrets
import shutil
import sqlite3
import os
import threading
import re
import zipfile

DATABASE_NAME = "users.db"  # go change in file_handler.py as well
CURRENT_FILE_DIRECTORY = os.path.dirname(os.path.abspath(__file__))  # home/sirat/Code/ADR-Prediction/WebBackend/
TYPING_TOOLS = ['hla_la', 'optitype', 'hisat_genotype', 'snp_bridge']
ASSOCIATED_FILETYPES = ['.bam', '.bam.bai', '.fq', '.sam']
SINGLE_USE_TOKENS = {}  # token used for downloading user data from server
UPLOAD_STREAMS = {}  # key: str(user.user_id) + sequence_label

app = Flask(__name__)
CORS(app)

# Increase file upload size limit
app.config['MAX_CONTENT_LENGTH'] = 300 * 1024 * 1024  # 300 MB


def remove_token(token):
    sleep(600)
    try:
        del SINGLE_USE_TOKENS[token]
    except KeyError:
        pass


def check_valid_label(name):
    max_size = 100
    return bool(re.fullmatch(r"[^\x00-\x1F\x7F]+", name)) and len(name) <= max_size


@app.route('/api/signup', methods=['POST'])
def api_signup_handler():
    data = request.json  # name email password
    log("Received signup data:", logtitle="POST REQUEST", var=data, color='yellow')

    user, return_message = User.preliminary_check(data)
    if return_message: return return_message, 400
    if user.add_user_to_db():
        return jsonify(UserSchema().dump(user, include=["user_id", "name", "email"]))

    return "Error during signup", 400


@app.route('/api/login', methods=['POST'])
def api_login_handler():
    data = request.json  # email password
    log("Received login data:", logtitle="POST REQUEST", var=data, color='yellow')

    user, return_message = User.preliminary_check(data)
    if return_message: return return_message, 400
    if not user.fetch_with_email_password(): return "Email or password is incorrect", 200

    user.create_session()
    return jsonify(UserSchema().dump(user, include=['user_id', 'name', 'email', 'cookie']))


@app.route('/api/checkEmail', methods=['POST'])
def api_check_same_email():
    data = request.json
    log("Received check email data:", logtitle="POST REQUEST", var=data, color='yellow')
    return jsonify({"value": check_value_in_column_exists(data["email"], "email")})


@app.route('/api/getAccountInfo', methods=['POST'])
def api_get_account_info():
    data = request.json  # cookie
    user, return_message = User.preliminary_check(data, ["cookie"])
    if return_message: return return_message, 400

    log("Requesting account info : ", logtitle="POST REQUEST", var=data, color='yellow')

    if user.fetch_with_cookie():
        return jsonify(UserSchema().dump(user, include=["user_id", "name", "email"]))
    else:
        return "Invalid cookie", 400


@app.route('/api/changeAccountInfo', methods=['POST'])
def api_change_account_info():
    data = request.json  # cookie + other params that you want to change
    log("Changing account info : ", logtitle="POST REQUEST", var=data, color='yellow')

    old_user = UserSchema().load({"cookie": data["cookie"]})
    if old_user.cookie is None: return "No session cookie provided", 400

    if not old_user.fetch_with_cookie():
        return "Invalid cookie", 400

    new_user = UserSchema().load(data)
    new_user.user_id = old_user.user_id

    if new_user.change_user_info(): return "All provided values changed", 200
    return "Invalid Information", 400


@app.route('/api/getSequences', methods=['POST'])
def api_get_sequences():
    data = request.json  # cookie
    user, return_message = User.preliminary_check(data, ["cookie"])
    if return_message: return return_message, 400

    log("Requesting sequences list : ", logtitle="POST REQUEST", var=data, color='yellow')

    connection, cursor = connect()
    cursor.execute("SELECT id,label,upload_time,status,igv,hla_la FROM User_sequences WHERE user_id = ?",
                   (user.user_id,))
    return jsonify({"list": cursor.fetchall()})


# Small file uploads ( <10 MB) and Large file uploads can use the same function.
@app.route("/api/uploadFile", methods=["POST"])
def api_upload_file():
    if 'chunk' not in request.files: return 'No file chunk', 400
    chunk = request.files['chunk']
    sequence_label = chunk.filename  # Label
    chunk_index = int(request.form['chunk_index'])
    total_chunks = int(request.form['total_chunk'])
    is_paired = request.form["is_paired"].lower() == "true"  # Me and the bois hate FormData

    log("Requesting to upload file (chunk) : ", logtitle="POST REQUEST", var=chunk_index, color='yellow')

    # get account info
    data = {"cookie": request.headers.get('Authorization').split('Bearer ')[-1]}
    user = UserSchema().load(data)
    if user.cookie is None: return "No session cookie provided", 400
    if not user.fetch_with_cookie(): return "Invalid cookie", 400
    if not check_valid_label(
            sequence_label): return "Invalid label (Genuinely tho how did this happen? I gave you emojis and you treat me like this?)", 400

    upload_identifier = str(
        user.user_id) + sequence_label  # I'm starting to realize that this could theoretically be a problem
    if chunk_index == 0:  # Initialize a new upload stream
        connection, cursor = connect()
        cursor.execute("INSERT INTO User_sequences (user_id, label, is_paired) VALUES (?,?,?)",
                       (user.user_id, sequence_label, is_paired))
        cursor.execute("SELECT id FROM User_sequences WHERE user_id = ? AND label = ?", (user.user_id, sequence_label))
        sequence_id = cursor.fetchall()[-1][0]
        disconnect(connection, cursor)
        UPLOAD_STREAMS[upload_identifier] = [sequence_id,
                                             []]  # empty list is list of uploaded chunks (to prevent shenanigans)

    sequence_id = str(UPLOAD_STREAMS[upload_identifier][0])

    if chunk_index in UPLOAD_STREAMS[upload_identifier][1]:
        log("MAYDAY MAYDAY SOMEONE MIGHT BE DOING SHENANIGANS", logtitle="error", var=sequence_id, color="red")
        return "Currently processing upload with same label please choose another.", 400

    UPLOAD_STREAMS[upload_identifier][1].append(chunk_index)

    log(f"Uploading to sequence id :", var=sequence_id)
    if not os.path.exists(user.path):
        os.makedirs(user.path)
        log(f"Directory '{user.path}' created.")

    directory = os.path.join(user.path, sequence_id)
    if not os.path.exists(directory):
        os.makedirs(directory)
        log(f"Directory '{directory}' created.")

    directory = os.path.join(directory, "chunk")
    if not os.path.exists(directory):
        os.makedirs(directory)
        log(f"Directory '{directory}' created.")

    chunk_path = os.path.join(directory, f"chunk_{chunk_index}")
    chunk.save(chunk_path)

    if len(os.listdir(directory)) != total_chunks:
        return f'Chunk {chunk_index} received', 200

    if not is_paired:

        final_path = os.path.join(user.path, f"{sequence_id}/{sequence_id}.fq")
        with open(final_path, 'ab') as final_file:
            for i in range(total_chunks):
                path_to_connect = os.path.join(directory, f'chunk_{i}')
                with open(path_to_connect, 'rb') as chunk_file:
                    final_file.write(chunk_file.read())

    else:
        file_path = os.path.join(user.path, f"{sequence_id}/{sequence_id}.zip")
        with open(file_path, 'ab') as final_file:
            for i in range(total_chunks):
                path_to_connect = os.path.join(directory, f'chunk_{i}')
                with open(path_to_connect, 'rb') as chunk_file:
                    final_file.write(chunk_file.read())

        log("Preparing to unzip paired end reads", logtitle='INFO', var=sequence_id, color='blue')
        try:
            with zipfile.ZipFile(file_path) as file:
                zip_dir = os.path.join(user.path, f"{sequence_id}/zip")
                file.extractall(zip_dir)  # Extract all contents of file into zip
                for index, _filename in enumerate(os.listdir(zip_dir)):
                    # Rename the file to R1.fq, R2.fq
                    run(f"mv '{zip_dir}/{_filename}' '{user.path}/{sequence_id}/{sequence_id}.R{index + 1}.fq'",
                        directory)
        except zipfile.BadZipfile:
            log("Invalid zip file", logtitle="error", var=sequence_id, color="red")
            shutil.rmtree(os.path.join(user.path, sequence_id))
            return "Invalid zip file", 400

    log("Large file upload completed", logtitle='success', var=sequence_id, color='green')

    threading.Thread(target=transform_file_for_igv,
                     args=(f"uploads/{user.user_id}/{sequence_id}", sequence_id)).start()
    del UPLOAD_STREAMS[upload_identifier]
    return 'Large file uploaded successfully', 200


@app.route('/api/run/hla_la', methods=['POST'])
def api_run_hla_la():
    # get account info
    data = request.json  # cookie sequence_id
    user, return_message = User.preliminary_check(data, ["cookie", "sequence_id"])
    if return_message: return return_message, 400

    log("Requesting to run HLA-LA : ", logtitle="POST REQUEST", color='yellow')
    status = get_column(user.sequence_id, "hla_la")

    if status == "COMPLETED":
        return 'File already processed', 200
    elif status == "RUNNING":
        return 'File being processed', 200

    threading.Thread(target=run_hla_la,
                     args=(f"uploads/{user.user_id}/{user.sequence_id}", user.sequence_id)).start()
    return 'Running HLA-LA', 200


@app.route('/api/requestFile/igv', methods=['POST'])
def api_request_file_igv():
    # get account info
    data = request.json
    user, return_message = User.preliminary_check(data, ["cookie", "sequence_id"])
    if return_message: return return_message, 400

    log("Requesting file with label", var=user.sequence_id, logtitle='POST REQUEST', color='yellow')
    directory = f"{user.path}/{user.sequence_id}/igv"

    if not os.path.exists(directory):
        log("Path doesn't exist", logtitle="error", color='red')
        return "No file found", 400

    bam_directory = f"{directory}/{user.sequence_id}.bam"
    bam_bai_directory = f"{directory}/{user.sequence_id}.bam.bai"

    if not os.path.isfile(bam_directory) or not os.path.isfile(bam_bai_directory):
        return "File isn't ready for view yet", 504

    # creating tokens for the two files requested (.bam, .bam.bai)
    token_bam = secrets.token_urlsafe(16)
    token_bam_bai = secrets.token_urlsafe(16)
    SINGLE_USE_TOKENS[token_bam] = bam_directory
    SINGLE_USE_TOKENS[token_bam_bai] = bam_bai_directory

    # remove the tokens in 10 minutes parallel
    threading.Thread(target=remove_token, args=(token_bam,)).start()
    threading.Thread(target=remove_token, args=(token_bam_bai,)).start()

    return jsonify({"token_bam": token_bam,
                    "token_bam_bai": token_bam_bai,
                    "range": get_range_of_bam_file(bam_directory)})


@app.route('/api/getResults/hla_la', methods=['POST'])
def api_get_typing_results():
    # get account info
    data = request.json
    user, return_message = User.preliminary_check(data, ["cookie", "sequence_id"])
    if return_message: return return_message, 400

    log("Requesting typing results with label :", var=user.sequence_id, logtitle='POST REQUEST', color='yellow')

    if get_column(user.sequence_id, "hla_la") != "COMPLETED":
        return jsonify([]), 200

    alleles = format_hla_la(f"{user.path}/{user.sequence_id}/hla_la/out/hla/R1_bestguess_G.txt")
    return_list = []
    for allele in alleles:
        return_list.extend(get_adr(allele))

    return jsonify(return_list)


@app.route('/api/deleteSequence', methods=['POST'])
def api_delete_sequence():
    # get account info
    data = request.json
    user, return_message = User.preliminary_check(data, ["cookie", "sequence_id"])
    if return_message: return return_message, 400

    log(f"Receiving delete sequence request with label :", var=user.sequence_id, logtitle="post request",
        color='yellow')

    directory = f"{user.path}/{user.sequence_id}"
    if os.path.exists(directory):
        shutil.rmtree(directory)

    log("Removed directory :", var=directory, logtitle="DELETE", color="red")

    connection, cursor = connect()
    cursor.execute("DELETE FROM User_sequences WHERE user_id=? AND id=?", (user.user_id, user.sequence_id))
    disconnect(connection, cursor)
    return "Successfully removed file", 200


@app.route("/api/getSequenceLog", methods=["POST"])
def api_get_sequence_log():  # I tried my best to use web sockets, and it didn't work. I give up man.
    # get account info
    data = request.json
    user, return_message = User.preliminary_check(data, ["cookie", "sequence_id"])
    if return_message: return return_message, 400

    log(f"Receiving request for logs", var=user.sequence_id, logtitle="post request", color='yellow')

    file_path = f"{user.path}/{user.sequence_id}/logs"

    if not os.path.exists(file_path):
        return jsonify({"logs": ""}), 200
    with open(file_path, 'r') as file:
        return jsonify({"logs": file.read()}), 200


@app.route("/data/sunburst", methods=["POST"])
def get_sunburst_data_using_filter():
    data = request.json
    countries = data["country_filter"]
    minP = data["minP_filter"]
    log('Receiving sunburst request with countries', var=countries, logtitle='post request', color='yellow')
    return jsonify(SunburstChartSetup.get_data(countries, minP))


@app.route("/statistic")
def get_statistic():
    connection, cursor = connect()
    cursor.execute("SELECT MAX(id) from User_sequences;")
    files_uploaded = cursor.fetchone()[0]
    disconnect(connection, cursor)
    return jsonify({"files_uploaded": files_uploaded}), 200


@app.route('/download/<token>')
def download_data_from_token(token):
    log(f"Receiving download request with token {token}", logtitle="download request", color='yellow')
    try:
        file_dir = SINGLE_USE_TOKENS[token]
        return send_file(file_dir)
    except KeyError:
        return "Invalid token", 400


@app.route('/reference/chr6.fa.fai')
def file_chr6_fa_fai():
    # Path to your sequence file
    sequence_file = 'references/chr6.fa.fai'

    # Use Flask's send_file function to send the file as a response
    return send_file(sequence_file)


@app.route('/reference/chr6.fa')
def file_chr6_fa():
    sequence_file = 'references/chr6.fa'

    return send_file(sequence_file)


@app.route('/reference/example_files/paired_end_example.zip')
def paired_end_example_zip():
    sequence_file = 'references/example_files/paired_end_example.zip'
    return send_file(sequence_file)


@app.route('/reference/example_files/single_end_example.zip')
def single_end_example_zip():
    sequence_file = 'references/example_files/single_end_example.zip'
    return send_file(sequence_file)


@app.route('/')
def hello():
    return "<b>Just checking if the backend works or not lol nice port number am I right :)</b>"


if __name__ == '__main__':
    app.run(debug=True, port=7000, host="0.0.0.0")
