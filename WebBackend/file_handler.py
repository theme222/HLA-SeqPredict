from TechnicalToolsV2 import log
import pandas as pd
import os
import sqlite3
import pysam

hla_la_currently_running = False


def set_available():
    global hla_la_currently_running
    hla_la_currently_running = True


def is_available():
    return not hla_la_currently_running


def connect(database_name="users.db"):
    log(f"Connecting to database {database_name}")
    connection = sqlite3.connect("databases/" + database_name)
    return connection, connection.cursor()


def disconnect(connection, cursor):
    log("Disconnecting from database")
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
    a = cursor.fetchone()[0]
    disconnect(connection, cursor)
    return a


def run(command):
    log("Running command :", var=command, logtitle='run', color='cyan')
    os.system(command)


def transform_file_for_igv(directory, filename, user_id):
    set_status_code(1, "igv", filename, user_id)
    os.makedirs(f"/app/{directory}/igv")
    run(f"cp /app/{directory}/{filename}.fq /app/{directory}/igv/{filename}.fq")
    directory += "/igv"
    run(f"singularity exec -B /app:/tmp /app/singularity/bwa_0_7_18 bwa mem /tmp/references/chr6.fa '/tmp/{directory}/{filename}.fq' > '/app/{directory}/{filename}.sam'")
    run(f"singularity exec -B /app:/tmp /app/singularity/samtools_1_20 samtools sort '/tmp/{directory}/{filename}.sam' > '/app/{directory}/{filename}.bam'")
    run(f"singularity exec -B /app:/tmp /app/singularity/samtools_1_20 samtools index '/tmp/{directory}/{filename}.bam'")
    set_status_code(2, "igv", filename, user_id)


def run_hla_la(directory, filename, user_id):
    global hla_la_currently_running
    hla_la_currently_running = True
    try:
        set_status_code(1, "hla_la", filename, user_id)
        os.makedirs(f"/app/{directory}/hla_la")
        run(f"cp /app/{directory}/{filename}.fq /app/{directory}/hla_la/{filename}.fq")
        directory += "/hla_la"
        run(f'singularity exec -B /app:/tmp /app/singularity/bwa_0_7_18 bwa mem /tmp/references/hg19.fasta /tmp/{directory}/{filename}.fq > /app/{directory}/{filename}.sam')
        run(f'singularity exec -B /app:/tmp /app/singularity/samtools_1_20 samtools sort /tmp/{directory}/{filename}.sam > /app/{directory}/{filename}.bam')
        run(f'singularity exec -B /app:/tmp /app/singularity/samtools_1_20 samtools index /tmp/{directory}/{filename}.bam')
        run(f"singularity exec -B /app:/tmp /app/singularity/hlala104_3.sif /usr/local/bin/HLA-LA/src/HLA-LA.pl --BAM /tmp/{directory}/{filename}.bam --workingDir /tmp/{directory} --customGraphDir /tmp/references/graphs --graph PRG_MHC_GRCh38_withIMGT --sampleID out --maxThreads 12 --longReads pacbio")
        hla_la_currently_running = False
        set_status_code(2, "hla_la", filename, user_id)
    except Exception as e:
        hla_la_currently_running = False
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
    output = ["HLA-" + i for i in list(df['Allele'])]  # add the keyword HLA-
    for index, allele in enumerate(output):
        output[index] = allele[0:allele.find("*") + 6]  # return only up to the "*" + 5 char

    return output


def get_adr(allele: str):  # return [ { allele: str, drug: str, odds: str, adr: str, link: str } , ... ]
    """
    :param allele:
    :return:
    """

    # formatting the standard input example: HLA-A*01:02
    gene = allele.split('*')[0]
    allelegroup = allele.split('*')[1].split(':')[0]
    protein = allele.split('*')[1].split(':')[0]

    connection, cursor = connect(database_name='hla_adr.db')

    # Get allele id
    cursor.execute("SELECT allele_id FROM alleles WHERE gene=? AND allelegroup=? AND (protein=? OR protein IS NULL)",
                   (gene, allelegroup, protein))
    allele_id = cursor.fetchone()
    if allele_id is None:
        return []
    allele_id = allele_id[0]

    cursor.execute("SELECT drug_id,adr_id,odds_exposed,link FROM HLA_ADR WHERE allele_id=(?)",
                   (allele_id,))

    info = cursor.fetchall()
    return_list = []

    for hla_adr in info:
        stuff = {}
        cursor.execute("SELECT name FROM Drugs WHERE drug_id=(?)",
                       (hla_adr[0],))
        stuff["drug"] = cursor.fetchone()[0]
        cursor.execute("SELECT name FROM ADR WHERE adr_id=(?)",
                       (hla_adr[1],))
        stuff["adr"] = cursor.fetchone()[0]
        stuff["allele"] = allele
        stuff["odds"] = str(hla_adr[2])
        stuff["link"] = str(hla_adr[3])
        return_list.append(stuff)

    disconnect(connection,cursor)
    return return_list


if __name__ == '__main__':
    a = []
    for i in format_hla_la("uploads/16/example.txt"):
        a.extend(get_adr(i))
    print(a)
