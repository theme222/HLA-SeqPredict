import os
import sqlite3
import time
import pandas as pd
from TechnicalToolsV2 import log, time_convert

filename = "database-setup.sql"
database_name = "test-database.db"

# Delete existing database file
if os.path.exists(database_name):
    os.remove(database_name)

connection = sqlite3.connect(database_name)
cursor = connection.cursor()

csv_name = "HLA-ADR.csv"
df = pd.read_csv(csv_name)


def run_sql_script(sql_file_name):
    with open(sql_file_name, 'r') as sql_file:
        a = sql_file.read()

    cursor.executescript(a)


def get_value(column, table_name, find_dict):
    """
    The longest one liner in existence but it doesn't work tho lol (too difficult to debug)
    def get_value(column, table_name, find_dict): return cursor.execute(f"SELECT ? FROM ? WHERE " + " ".join([f"{key}={value} AND" for key, value in find_dict.items()]) + " TRUE;",(column,table_name)).fecthone()[0]
    """
    # use the keys and values in find_dict to generate a where clause
    find_string = ""
    for key, value in find_dict.items():
        if value is not None:
            find_string += f'{key}="{value}" AND '
        else:
            find_string += f"{key} IS NULL AND "

    # Find the item
    command = f"SELECT {column} FROM {table_name} WHERE " + find_string + "TRUE;"
    log(command, logtitle="running command", color='purple')
    cursor.execute(command)
    result = cursor.fetchone()
    return result[0] if result else None


def add_values():
    # Adding into allele table
    for i, e in enumerate(set(df["hla_adr_allele"])):
        gene = e.split("*")[0]
        group = e.split(":")[0].split("*")[1]
        try:
            protein = e.split(":")[1]
        except IndexError:
            protein = None
        print(i, gene, group, protein)
        cursor.execute(
            "INSERT INTO Alleles (allele_id, gene, allelegroup, protein) VALUES (?,?,?,?);",
            (i, gene, group, protein)
        )

    # Adding into drugs table
    for i, e in enumerate(set(df["hla_adr_drug_id"])):
        print(i, e)
        cursor.execute(f'INSERT INTO Drugs (drug_id, name) VALUES (?,?)', (i, e))

    # Adding into ADR table
    for i, e in enumerate(set(df["hla_adr_adr"])):
        print(i, e)
        cursor.execute(f"INSERT INTO ADR (adr_id, name) VALUES (?,?)", (i, e))

    # Adding into connections table
    alleles = list(df["hla_adr_allele"])
    drugs = list(df["hla_adr_drug_id"])
    adr = list(df["hla_adr_adr"])
    for i in range(len(adr)):

        # find the allele_id
        allele = alleles[i]
        gene = allele.split("*")[0]
        group = allele.split(":")[0].split("*")[1]
        try:
            protein = allele.split(":")[1]
        except IndexError:
            protein = None
        allele_id = get_value("allele_id", "Alleles", {"gene": gene, "allelegroup": group, "protein": protein})

        # find the drug_id
        drug_id = get_value("drug_id", "Drugs", {"name": drugs[i]})

        # find the adr_id
        adr_id = get_value("adr_id", "ADR", {"name": adr[i]})

        cursor.execute("INSERT INTO HLA_ADR (id,allele_id,drug_id,adr_id) VALUES (?,?,?,?)",
                       (i, allele_id, drug_id, adr_id))


def main():
    log("Running script ...", logtitle="START", color='cyan')
    start_time = time.time()
    run_sql_script(filename)

    add_values()

    connection.commit()
    cursor.close()
    connection.close()
    log(f"Total time elapsed : {time_convert(time.time() - start_time)}", logtitle="END", color="cyan")


if __name__ == '__main__':
    main()
