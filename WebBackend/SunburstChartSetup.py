import json
import time
from TechnicalToolsV2 import time_convert, log
from file_handler import CURRENT_FILE_DIRECTORY

'''
origin : item
child : list
key : string
value : string (maybe)
'''

file_path = f'{CURRENT_FILE_DIRECTORY}/databases/sunburst.json'
iteration = ["hla_adr_disease_id", "hla_adr_drug_id", "hla_adr_adr", "hla_adr_gene_id", "hla_adr_allele"]

with open(file_path, 'r') as file:
    data = json.load(file)  # Load JSON file into a Python dictionary


def find_index(child, value):
    for i, e in enumerate(child):
        if e["name"] == value:
            return i
    return -1


def child_recursion(origin, child, cycle, parent=None):  # Can't believe I am using recursion by my own free will
    if parent is None: parent = []

    value = origin[iteration[cycle]]
    if cycle == 4:
        if find_index(child, value) == -1:
            child.append({"name": value, "value": origin["GeTH_AlleleCount"], "parent": parent})
    else:
        if find_index(child, value) == -1:
            child.append({"name": value, "children": [], "parent": parent})
        index = find_index(child, value)
        new_parent = parent.copy()
        new_parent.append(value)
        child_recursion(origin, child[index]["children"], cycle + 1, new_parent)


def get_data(country, minP):
    final_list = []
    country = set(country)
    for item in data:
        if minP and item["hla_adr_minP"] == "FALSE":
            continue
        if item["hla_adr_country_id"] not in country:
            continue
        child_recursion(item, final_list, 0)

    return final_list


if __name__ == '__main__':
    start_time = time.time()
    stuff = get_data(["Thailand", "Canada"], False)
    with open("output.json",'w') as file:
        json.dump(stuff, file, indent=2)
    end_time = time.time()
    log("Completed in", var=time_convert(end_time - start_time))
