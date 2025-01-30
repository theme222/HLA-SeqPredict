import json
import time
from TechnicalToolsV2 import log, time_convert

filename = "sunburst.json"
outname = "sunburst.table.json"

"""
disease
drug
allele
gene_loci
adr
ethnicity
country
pvalue_exposed
pvalue_population
odds_exposed
pubmed
allele_count
allele_freq
minP
"""

with open(filename) as file:
    data = json.load(file)

def move_decimal_right(num):  # I am not sorry
    num_str = str(num)
    stuff = num_str.split(".")
    int_part = stuff[0]
    float_part = stuff[1]
    int_part += float_part[:2]
    float_part = float_part[2:]
    return float(int_part + '.' + float_part)

def get_value(key,item):
    if key == "hla_adr_minP": return item[key] == "TRUE"
    if key in item: return item[key]
    else: return ""

def main():
    output = []
    for item in data:
        if item["GeTH_AlleleCount"] == 0:
            continue
        obj = {
            "disease": get_value("hla_adr_disease_id", item),
            "drug": get_value("hla_adr_drug_id",item),
            "allele": get_value("hla_adr_allele", item),
            "gene_loci": get_value("hla_adr_gene_id", item),
            "adr": get_value("hla_adr_adr", item),
            "ethnicity": get_value("hla_adr_cohort_ethnicity_id", item),
            "country": get_value("hla_adr_country_id", item),
            "pvalue_exposed": get_value("hla_adr_pvalue_exposed", item),
            "pvalue_population": get_value("hla_adr_pvalue_population", item),
            "odds_exposed": get_value("hla_adr_odds_exposed", item),
            "pubmed": get_value("hla_adr_pubmed_link", item),
            "allele_count": get_value("GeTH_AlleleCount", item),
            "allele_freq": move_decimal_right(get_value("GeTH_AlleleFreq", item)),
            "minP" : get_value("hla_adr_minP", item)
        }
        output.append(obj)



    with open(outname, 'w') as file:
        json.dump(output, file, indent=2)

if __name__ == '__main__':
    start_time = time.time()
    main()
    log("Completed in", var=time_convert(time.time() - start_time))