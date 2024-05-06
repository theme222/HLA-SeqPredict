from TechnicalToolsV2 import log
import pandas as pd
import os
import sqlite3
import pysam


def connect(database_name="users.db"):
    log(f"Connecting to database {database_name}")
    connection = sqlite3.connect(database_name)
    return connection, connection.cursor()


def disconnect(connection, cursor):
    log("Disconnecting from database")
    connection.commit()
    cursor.close()
    connection.close()


def run(command):
    log("Running command :", var=command, logtitle='run', color='cyan')
    os.system(command)


def transform_file_for_igv(directory, filename, user_id):
    run(f"bwa mem references/chr6.fa '{directory}/{filename}.fq' > '{directory}/{filename}.sam'")
    run(f"samtools sort '{directory}/{filename}.sam' > '{directory}/{filename}.bam'")
    run(f"samtools index '{directory}/{filename}.bam'")

    connection, cursor = connect()
    cursor.execute("UPDATE User_sequences SET igv = 1 WHERE user_id = ? AND label = ?", (user_id, filename))
    disconnect(connection, cursor)


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


def format_hla_la(user_id, label, path_to_output):
    connection, cursor = connect()

    # Check if hlala is done processing
    cursor.execute("SELECT hlala FROM User_sequences WHERE user_id=? AND label=?", (user_id, label))
    if cursor.fetchall()[0] == 0:
        return None

    disconnect(connection, cursor)

    # format the output
    df = pd.read_csv(path_to_output, sep='\t')
    output = ["HLA-" + i for i in list(df['Allele'])]  # add the keyword HLA-
    for index, allele in enumerate(output):
        output[index] = allele[0:allele.find("*") + 6]  # return only up to the "*" + 5 char

        gene = output[index].split('*')[0]
        allelegroup = output[index].split('*')[1].split(':')[0]
        protein = output[index].split('*')[1].split(':')[0]
        output[index] = {"gene": gene, "allelegroup": allelegroup, "protein": protein}

    return output

def get_adr(formatted_output: dict):
    """
    :param formatted_output:
    :return:
    """
    connection, cursor = connect(database_name='hla_adr.db')
    cursor.execute("SELECT name FROM drugs WHERE drug_id="
                   "(SELECT drug_id FROM HLA_ADR WHERE allele_id=("
                   "SELECT allele_id FROM alleles WHERE gene=? AND allelegroup=? AND (protein=? OR protein IS NULL)))",
                   (formatted_output['gene'], formatted_output['allelegroup'], formatted_output['protein']))
    a = cursor.fetchall()
    disconnect(connection,cursor)
    return a


if __name__ == '__main__':
    for i in format_hla_la(16, "file1", "uploads/16/example.txt"):
        print(get_adr(i))
