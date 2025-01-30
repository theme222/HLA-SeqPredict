import json
from TechnicalToolsV2 import log
from flask import jsonify
from marshmallow import Schema, fields, post_load, EXCLUDE
from time import time, sleep
import secrets
import pandas as pd
import os
import sqlite3
import threading
import pysam

CURRENT_FILE_DIRECTORY = os.path.dirname(os.path.abspath(__file__))
HLA_ADR_PATH = 'databases/hla_adr.json'
SESSION_DURATION = 1000000  # Seconds
DATABASE_TIMEOUT = 0.5  # seconds


class User:
    def __init__(self, user_id, name, email, password, cookie):
        # Use "UserSchema().load(data)" to load this user | data : dict(str: any)
        self.user_id = user_id
        self.name = name
        self.email = email
        self.password = password
        self.cookie = cookie

    def add_user_to_db(self):
        log("Signing up user : ", logtitle="CLASS METHOD", var=self.name, color='black')

        if self.name is None or self.email is None or self.password is None:
            log("Required information doesn't exist", logtitle="error", var=self.email, color="red")
            return False
        if check_value_in_column_exists(self.email, "email"):
            log("User already exists", logtitle="error", var=self.email, color='red')
            return False

        connection, cursor = connect()
        cursor.execute(
            "INSERT INTO Users (name, email, password_hash) VALUES (?,?,?);",
            (self.name, self.email, self.password)
        )
        disconnect(connection, cursor)
        self.fetch_with_email_password()
        return True

    def create_session(self):
        self.cookie = generate_session_cookie()

        connection, cursor = connect()
        cursor.execute(f"INSERT INTO Sessions (user_id,token,start_time,duration) VALUES (?,?,?,?)",
                       (self.user_id, self.cookie, int(time()), SESSION_DURATION))
        disconnect(connection, cursor)

        log("Created cookie : ", var=self.cookie, color="green", logtitle="COOKIE MONSTER")

    def fetch_with_email_password(self):
        log("Fetching email : ", logtitle="CLASS METHOD", var=self.email, color='black')
        connection, cursor = connect()
        cursor.execute("SELECT id,name,email,password_hash From Users Where email = ? AND password_hash = ?",
                       (self.email, self.password))
        info = cursor.fetchone()
        disconnect(connection, cursor)

        if info:
            self.user_id, self.name, self.email, self.password = info
            return True
        else:
            log("Invalid email or password", logtitle='error', color='red')
            return False

    def fetch_with_cookie(self, get_password=False, change_values=True):
        log("Fetching cookie : ", logtitle="CLASS METHOD", var=self.cookie, color='black')
        connection, cursor = connect()
        cursor.execute(
            "SELECT id,name,email,password_hash FROM Users WHERE id = (SELECT user_id FROM Sessions WHERE token = ?)",
            (self.cookie,))
        info = cursor.fetchone()
        disconnect(connection, cursor)

        if info:
            if change_values: self.user_id, self.name, self.email, self.password = info
            if not get_password: self.password = None  # Probably a vulnerability here but I don't really care
            return True
        else:
            log("Invalid cookie", logtitle='error', color='red')
            return False

    def change_user_info(self):
        """
        Runs on new User object with correct user_id and any new information wanting to be changed.
        Any property left blank (None) will not be changed.
        This will change the info inside the database.
        :return bool:
        """
        log("Changing info of user_id : ", logtitle="CLASS METHOD", var=self.user_id, color='black')

        # Checks for valid user_id
        if self.user_id is None:
            log("No user_id provided", logtitle="error", color="red")
            return False
        if not check_value_in_column_exists(self.user_id, "id"):
            log("User not found", logtitle="error", color="red")
            return False

        # Checks for valid email
        if check_value_in_column_exists(self.email, "email"):
            log("Email already in use", logtitle="error", color='red')
            return False

        connection, cursor = connect()
        for key, value in UserSchema().dump(self, include=["name", "email", "password"]).items():
            if value is not None:
                if key == "password": key += "_hash"  # I wanna cry
                cursor.execute(f"UPDATE Users SET {key}=? WHERE id=?", (value, self.user_id))

        disconnect(connection, cursor)

        return True

    @property
    def path(self):
        return f"{CURRENT_FILE_DIRECTORY}/uploads/{self.user_id}"

    def __bool__(self):  # Use this function to check if user exists
        return self.user_id is not None


class UserSchema(Schema):
    user_id = fields.Int(missing=None)
    name = fields.Str(missing=None)
    email = fields.Email(missing=None)
    password = fields.Str(missing=None)
    cookie = fields.Str(missing=None)

    class Meta:
        unknown = EXCLUDE

    @post_load
    def make_user(self, data, **kwargs):
        return User(**data)

    def dump(self, obj,
             include=None):  # Untested chatgpt code example : schema.dump(userObject, include=['user_id','password']
        if include:
            data = {key: getattr(obj, key) for key in include if hasattr(obj, key)}
            return super().dump(data)
        return super().dump(obj)


def generate_session_cookie():
    return f"{secrets.token_hex(16)}__{int(time())}"


def check_value_in_column_exists(value, column_name):
    connection, cursor = connect()
    cursor.execute(f"SELECT EXISTS(SELECT 1 FROM Users WHERE {column_name} = ?) ", (value,))
    exists = cursor.fetchone()[0]
    disconnect(connection, cursor)
    return exists == 1


def connect(database_name="users.db"):
    log(f"Connecting to database {database_name}", logtitle="database", color='black')
    connection = sqlite3.connect("databases/" + database_name)
    cursor = connection.cursor()
    threading.Thread(target=disconnect, args=(connection, cursor, True)).start()
    return connection, cursor


def disconnect(connection, cursor, _internalCall=False):
    if _internalCall:  # This is for when an error happens while executing sql commands
        sleep(DATABASE_TIMEOUT)
        try:
            connection.commit()
            cursor.close()
            connection.close()
            log("Automatically disconnecting from database due to timeout", logtitle="timeout", color='red',
                var=DATABASE_TIMEOUT)
        except sqlite3.ProgrammingError:
            pass
        finally:
            return

    log("Disconnecting from database", logtitle="database", color='black')
    connection.commit()
    cursor.close()
    connection.close()


def set_status_code(value_to_set, column_name, label, user_id):
    connection, cursor = connect()
    cursor.execute(f"UPDATE User_sequences SET {column_name} = ? WHERE label = ? AND user_id = ?",
                   (value_to_set, label, user_id))
    disconnect(connection, cursor)


def get_status_code(label, column_name, user_id):
    connection, cursor = connect()
    cursor.execute(f"SELECT {column_name} FROM User_sequences WHERE label = ? AND user_id = ?", (label, user_id))
    _a = cursor.fetchone()[0]
    disconnect(connection, cursor)
    return _a


def run(command):
    log("Running command :", var=command, logtitle='run', color='cyan')
    os.system(command)


def transform_file_for_igv(directory, filename, user_id):
    set_status_code(1, "igv", filename, user_id)
    os.makedirs(f"/app/{directory}/igv")

    # TODO: Make the checking for if it is paired end related to the database
    if os.path.exists(f"{directory}/zip"):
        run(f"cp /app/{directory}/{filename}.R1.fq /app/{directory}/igv/{filename}.R1.fq")
        run(f"cp /app/{directory}/{filename}.R2.fq /app/{directory}/igv/{filename}.R2.fq")
        directory += "/igv"
        run(f"singularity exec -B /app:/tmp /app/singularity/bwa_0_7_18 bwa mem /tmp/references/chr6.fa '/tmp/{directory}/{filename}.R1.fq' '/tmp/{directory}/{filename}.R2.fq' > '/app/{directory}/{filename}.sam'")
    else:
        run(f"cp /app/{directory}/{filename}.fq /app/{directory}/igv/{filename}.fq")
        directory += "/igv"
        run(f"singularity exec -B /app:/tmp /app/singularity/bwa_0_7_18 bwa mem /tmp/references/chr6.fa '/tmp/{directory}/{filename}.fq' > '/app/{directory}/{filename}.sam'")

    run(f"singularity exec -B /app:/tmp /app/singularity/samtools_1_20 samtools sort '/tmp/{directory}/{filename}.sam' > '/app/{directory}/{filename}.bam'")
    run(f"singularity exec -B /app:/tmp /app/singularity/samtools_1_20 samtools index '/tmp/{directory}/{filename}.bam'")
    set_status_code(2, "igv", filename, user_id)


def run_hla_la(directory, filename, user_id):
    try:
        set_status_code(1, "hla_la", filename, user_id)
        os.makedirs(f"/app/{directory}/hla_la")

        # TODO: Make the checking for if it is paired end related to the database
        if os.path.exists(f"{directory}/zip"):
            run(f"cp /app/{directory}/{filename}.R1.fq /app/{directory}/hla_la/{filename}.R1.fq")
            run(f"cp /app/{directory}/{filename}.R2.fq /app/{directory}/hla_la/{filename}.R2.fq")
            directory += "/hla_la"
            run(f'singularity exec -B /app:/tmp /app/singularity/bwa_0_7_18 bwa mem /tmp/references/hg19.fasta /tmp/{directory}/{filename}.R1.fq /tmp/{directory}/{filename}.R2.fq > /app/{directory}/{filename}.sam')
            run(f'singularity exec -B /app:/tmp /app/singularity/samtools_1_20 samtools sort /tmp/{directory}/{filename}.sam > /app/{directory}/{filename}.bam')
            run(f'singularity exec -B /app:/tmp /app/singularity/samtools_1_20 samtools index /tmp/{directory}/{filename}.bam')
            run(f"singularity exec -B /app:/tmp /app/singularity/hlala104_3.sif /usr/local/bin/HLA-LA/src/HLA-LA.pl --BAM /tmp/{directory}/{filename}.bam --workingDir /tmp/{directory} --customGraphDir /tmp/references/graphs --graph PRG_MHC_GRCh38_withIMGT --sampleID out --maxThreads 32")
        else:
            run(f"cp /app/{directory}/{filename}.fq /app/{directory}/hla_la/{filename}.fq")
            directory += "/hla_la"
            run(f'singularity exec -B /app:/tmp /app/singularity/bwa_0_7_18 bwa mem /tmp/references/hg19.fasta /tmp/{directory}/{filename}.fq > /app/{directory}/{filename}.sam')
            run(f'singularity exec -B /app:/tmp /app/singularity/samtools_1_20 samtools sort /tmp/{directory}/{filename}.sam > /app/{directory}/{filename}.bam')
            run(f'singularity exec -B /app:/tmp /app/singularity/samtools_1_20 samtools index /tmp/{directory}/{filename}.bam')
            run(f"singularity exec -B /app:/tmp /app/singularity/hlala104_3.sif /usr/local/bin/HLA-LA/src/HLA-LA.pl --BAM /tmp/{directory}/{filename}.bam --workingDir /tmp/{directory} --customGraphDir /tmp/references/graphs --graph PRG_MHC_GRCh38_withIMGT --sampleID out --maxThreads 32 --longReads ont2d")

        set_status_code(2, "hla_la", filename, user_id)
    except Exception as e:
        set_status_code(0, "hla_la", filename, user_id)
        raise e


def get_range_of_bam_file(bam_file):
    # Open the BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Iterate through alignments
        first_position = None
        last_position = None
        for read in bam:
            if first_position is None:
                first_position = read.reference_start
            last_position = read.reference_end
        return f"{first_position}-{last_position}"


def format_hla_la(path_to_output):
    """
    if get_status_code(label, "hla_la", user_id) != 2:  # It isn't completed
        return None
    """
    # format the output
    df = pd.read_csv(path_to_output, sep='\t')
    print(df.keys())
    output = ["HLA-" + i for i in list(df['Allele']) if
              i[0] in ["A", "B", 'C']]  # add the keyword HLA- and only focus A, B, C
    for index, allele in enumerate(output):
        output[index] = allele[0:allele.find("*") + 6]  # return only up to the "*" + 5 char

    return output


def get_adr(allele: str):  # return [ { allele: str, drug: str, odds_exposed: str, adr: str, pubmed: str } , ... ]
    """
    :param allele:
    :return:
    """

    # formatting the standard input example: HLA-A*01:02
    gene = allele.split('*')[0]
    allelegroup = allele.split('*')[1].split(':')[0]
    protein = allele.split('*')[1].split(':')[0]

    with open(HLA_ADR_PATH, "r") as file:
        file = json.load(file)

    return_list = []

    for item in file:
        if item["allele"] == allele:
            return_list.append(item)

    return return_list


if __name__ == '__main__':
    a = []
    for i in format_hla_la("uploads/16/example.txt"):
        a.extend(get_adr(i))
    print(a)
