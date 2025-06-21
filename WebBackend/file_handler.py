import json
from TechnicalToolsV2 import log
from marshmallow import Schema, fields, post_load, EXCLUDE, ValidationError
from time import time, sleep
import secrets
import pandas as pd
import os
import sqlite3
import threading
import pysam
import subprocess
import sys

CURRENT_FILE_DIRECTORY = os.path.dirname(os.path.abspath(__file__))
HLA_ADR_PATH = 'databases/hla_adr.json'
SESSION_DURATION = 1000000  # Seconds
DATABASE_TIMEOUT = 0.5  # seconds


class User:
    def __init__(self, user_id: int, name: str, email: str, password: str, cookie: str, sequence_id: int):
        # Use "UserSchema().load(data)" to load this user | data : dict(str: any)
        self.user_id = user_id
        self.name = name
        self.email = email
        self.password = password
        self.cookie = cookie
        self.sequence_id = sequence_id

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

    def fetch_with_cookie(self, get_password=False, change_values=True) -> bool:
        log("Fetching cookie : ", logtitle="CLASS METHOD", var=self.cookie, color='black')

        connection, cursor = connect()

        cursor.execute("SELECT user_id, start_time, duration From Sessions WHERE token = ?", (self.cookie,))
        session_info = cursor.fetchone()

        if not session_info: return False
        user_id, start_time, duration = session_info
        if int(time()) > start_time + duration:
            cursor.execute("DELETE FROM Sessions WHERE token = ?", (self.cookie,))
            return False

        cursor.execute(
            "SELECT id,name,email,password_hash FROM Users WHERE id = ?",
            (user_id,))
        info = cursor.fetchone()

        disconnect(connection, cursor)

        if info:
            if change_values: self.user_id, self.name, self.email, self.password = info
            if not get_password: self.password = None  # Probably a vulnerability here but I don't really care
            return True
        else:
            log("User doesn't exist", logtitle='error', color='red')
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

    def validate_sequence_ownership(self) -> bool:
        """
        :return bool:
        Pass through the sequence id, and it will check if the user owns the sequence
        """
        connection, cursor = connect()
        cursor.execute(f"SELECT user_id FROM User_sequences WHERE id=?", (self.sequence_id,))
        sequence_owner_id = cursor.fetchone()
        disconnect(connection, cursor)
        return sequence_owner_id is None or self.user_id == sequence_owner_id[0]

    @staticmethod
    def preliminary_check(data, properties_to_check: list[str] = None):
        """
        :param data:
        :param properties_to_check: A list of user properties that are expected to exist and be validated. example: [user_id, name, email, cookie]
        :return:
        Does common checks while loading user data and returns a user object and an error code if it failed a check
        """
        if properties_to_check is None:
            properties_to_check = []

        try:
            user: User = UserSchema().load(data)
        except ValidationError:
            return UserSchema().load({}), "Provided values are in an incorrect format"

        if "cookie" in properties_to_check and user.cookie is None:
            return user, "No session cookie provided"
        if "cookie" in properties_to_check and not user.fetch_with_cookie():
            return user, "Invalid cookie"
        if "sequence_id" in properties_to_check and not user.validate_sequence_ownership():
            return user, "Sequence doesn't belong to user"

        return user, ""

    @property
    def path(self):
        return f"{CURRENT_FILE_DIRECTORY}/uploads/{self.user_id}"

    def __bool__(self):  # Use this function to check if user exists
        return self.user_id is not None


class UserSchema(Schema):
    user_id = fields.Int(load_default=None)
    name = fields.Str(load_default=None)
    email = fields.Email(load_default=None)
    password = fields.Str(load_default=None)
    cookie = fields.Str(load_default=None)
    sequence_id = fields.Int(load_default=None)

    class Meta:
        unknown = EXCLUDE

    @post_load
    def make_user(self, data, **kwargs):
        return User(**data)

    def dump(self, obj, include=None):  # Untested chatgpt code example : schema.dump(userObject, include=['user_id','password']
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
    # TODO: Consider the pros and cons of adding a context manager to this section instead of the connect disconnect bullshit
    log(f"Connecting to database {database_name}", logtitle="database", color='black')
    connection = sqlite3.connect("databases/" + database_name)
    cursor = connection.cursor()
    threading.Thread(target=disconnect, args=(connection, cursor, True)).start()  # This is the internal Call
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


def set_column(value_to_set, column_name, sequence_id, table="User_sequences"):
    connection, cursor = connect()
    cursor.execute(f"UPDATE {table} SET {column_name} = ? WHERE id = ?",
                   (value_to_set, sequence_id))
    disconnect(connection, cursor)


def get_column(sequence_id, column_name, table="User_sequences"):
    connection, cursor = connect()
    cursor.execute(f"SELECT {column_name} FROM {table} WHERE id = ?", (sequence_id,))
    value = cursor.fetchone()[0]
    disconnect(connection, cursor)
    return value


def run(command, log_directory=None):
    """
    :param command:
    :param log_directory:
    :return:

    The values being written into the file is as follows
    > command being executed
    ! error
    logsandmessages
    """
    log("Running command :", var=command, logtitle='run', color='cyan')

    if log_directory is None:
        os.system(command)
        return

    file_path = os.path.join(log_directory, "logs")

    file = open(file_path, 'a')
    try:
        file.write(f"> {command}\n")
        # Process has stdout and stderr which will buffer for each new line thus giving "real time" updates hopefully
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                                   bufsize=1)

        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            file.write(line)
            file.flush()

        for line in process.stderr:
            sys.stderr.write(line)
            sys.stderr.flush()
            file.write(line)
            file.flush()

        process.wait()

    except Exception as e:
        if file: file.write("! " + str(e))
        raise e

    finally:
        if file: file.close()


def transform_file_for_igv(directory, sequence_id):
    # directory : uploads/user_id/sequence_id

    set_column("RUNNING", "igv", sequence_id)
    set_column("RUNNING IGV", "status", sequence_id)
    os.makedirs(f"/app/{directory}/igv")

    sequence_directory = f"/app/{directory}"
    is_paired = get_column(sequence_id, "is_paired") == 1  # sqlite does not have a bool value

    try:
        if is_paired:
            run(f"cp /app/{directory}/{sequence_id}.R1.fq /app/{directory}/igv/{sequence_id}.R1.fq", sequence_directory)
            run(f"cp /app/{directory}/{sequence_id}.R2.fq /app/{directory}/igv/{sequence_id}.R2.fq", sequence_directory)
            directory += "/igv"
            run(f"singularity exec -B /app:/tmp /app/singularity/bwa_0_7_18 bwa mem /tmp/references/chr6.fa '/tmp/{directory}/{sequence_id}.R1.fq' '/tmp/{directory}/{sequence_id}.R2.fq' > '/app/{directory}/{sequence_id}.sam'",
                sequence_directory)
        else:
            run(f"cp /app/{directory}/{sequence_id}.fq /app/{directory}/igv/{sequence_id}.fq", sequence_directory)
            directory += "/igv"
            run(f"singularity exec -B /app:/tmp /app/singularity/bwa_0_7_18 bwa mem /tmp/references/chr6.fa '/tmp/{directory}/{sequence_id}.fq' > '/app/{directory}/{sequence_id}.sam'",
                sequence_directory)

        run(f"singularity exec -B /app:/tmp /app/singularity/samtools_1_20 samtools sort '/tmp/{directory}/{sequence_id}.sam' > '/app/{directory}/{sequence_id}.bam'",
            sequence_directory)
        run(f"singularity exec -B /app:/tmp /app/singularity/samtools_1_20 samtools index '/tmp/{directory}/{sequence_id}.bam'",
            sequence_directory)
        set_column("COMPLETED", "igv", sequence_id)
    except Exception as e:
        set_column("ERROR", "igv", sequence_id)
        raise e
    finally:
        set_column("IDLE", "status", sequence_id)


def run_hla_la(directory, sequence_id):
    # directory : uploads/user_id/sequence_id
    set_column("RUNNING", "hla_la", sequence_id)
    set_column("RUNNING HLA*LA", "status", sequence_id)

    sequence_directory = f"/app/{directory}"

    is_paired = get_column(sequence_id, "is_paired") == 1

    try:
        os.makedirs(f"/app/{directory}/hla_la")

        if is_paired:
            run(f"cp /app/{directory}/{sequence_id}.R1.fq /app/{directory}/hla_la/{sequence_id}.R1.fq",
                sequence_directory)
            run(f"cp /app/{directory}/{sequence_id}.R2.fq /app/{directory}/hla_la/{sequence_id}.R2.fq",
                sequence_directory)
            directory += "/hla_la"
            run(f'singularity exec -B /app:/tmp /app/singularity/bwa_0_7_18 bwa mem /tmp/references/hg19.fasta /tmp/{directory}/{sequence_id}.R1.fq /tmp/{directory}/{sequence_id}.R2.fq > /app/{directory}/{sequence_id}.sam',
                sequence_directory)
            run(f'singularity exec -B /app:/tmp /app/singularity/samtools_1_20 samtools sort /tmp/{directory}/{sequence_id}.sam > /app/{directory}/{sequence_id}.bam',
                sequence_directory)
            run(f'singularity exec -B /app:/tmp /app/singularity/samtools_1_20 samtools index /tmp/{directory}/{sequence_id}.bam',
                sequence_directory)
            run(f"singularity exec -B /app:/tmp /app/singularity/hlala104_3.sif /usr/local/bin/HLA-LA/src/HLA-LA.pl --BAM /tmp/{directory}/{sequence_id}.bam --workingDir /tmp/{directory} --customGraphDir /tmp/references/graphs --graph PRG_MHC_GRCh38_withIMGT --sampleID out --maxThreads 32",
                sequence_directory)
        else:
            run(f"cp /app/{directory}/{sequence_id}.fq /app/{directory}/hla_la/{sequence_id}.fq", sequence_directory)
            directory += "/hla_la"
            run(f'singularity exec -B /app:/tmp /app/singularity/bwa_0_7_18 bwa mem /tmp/references/hg19.fasta /tmp/{directory}/{sequence_id}.fq > /app/{directory}/{sequence_id}.sam',
                sequence_directory)
            run(f'singularity exec -B /app:/tmp /app/singularity/samtools_1_20 samtools sort /tmp/{directory}/{sequence_id}.sam > /app/{directory}/{sequence_id}.bam',
                sequence_directory)
            run(f'singularity exec -B /app:/tmp /app/singularity/samtools_1_20 samtools index /tmp/{directory}/{sequence_id}.bam',
                sequence_directory)
            run(f"singularity exec -B /app:/tmp /app/singularity/hlala104_3.sif /usr/local/bin/HLA-LA/src/HLA-LA.pl --BAM /tmp/{directory}/{sequence_id}.bam --workingDir /tmp/{directory} --customGraphDir /tmp/references/graphs --graph PRG_MHC_GRCh38_withIMGT --sampleID out --maxThreads 32 --longReads ont2d",
                sequence_directory)

        set_column("COMPLETED", "hla_la", sequence_id)
    except Exception as e:
        set_column("ERROR", "hla_la", sequence_id)
        raise e
    finally:
        set_column("IDLE", "status", sequence_id)


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
