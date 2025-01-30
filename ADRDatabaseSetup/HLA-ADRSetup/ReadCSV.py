import pandas as pd

# Read the file into a pandas DataFrame
df = pd.read_csv('HLA-ADR.csv')

# Now you can work with the DataFrame as needed
a = set(df["hla_adr_allele"])
for i,e in enumerate(a):
    gene = e.split("*")[0]
    group = e.split(":")[0].split("*")[1]
    try:
        protein = e.split(":")[1]
    except IndexError:
        protein = None
    print(i,gene, group, protein)
