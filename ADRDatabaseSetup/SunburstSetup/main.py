import json
import time
from TechnicalToolsV2 import time_convert, log

'''
origin : item
child : list
key : string
value : string (maybe)
'''


file_path = 'sunburst.json'
iteration = ["hla_adr_disease_id", "hla_adr_drug_id", "hla_adr_adr", "hla_adr_gene_id", "hla_adr_allele"]
colors = ["#DB6F57", "#ECE0C5", "#CF648C", "#5DA597", "#55547F"]
out_path = file_path + ".out.json"

with open(file_path, 'r') as file:
    data = json.load(file)  # Load JSON file into a Python dictionary


def find_index(child, value):
    for i, e in enumerate(child):
        if e["name"] == value:
            return i
    return -1



def thing_idk(origin, child, cycle, parent=None):  # Can't believe I am using recursion by my own free will
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
        thing_idk(origin, child[index]["children"], cycle + 1, new_parent)


def main():
    disease = set()
    _disease = "hla_adr_disease_id"
    drug = set()
    _drug = "hla_adr_drug_id"
    adr = set()
    _adr = "hla_adr_adr"
    group = set()
    _group = "hla_adr_gene_id"
    allele = set()
    _allele = "hla_adr_allele"

    for item in data:
        disease.add(item[_disease])
        drug.add(item[_drug])
        adr.add(item[_adr])
        group.add(item[_group])
        allele.add(item[_allele])

    final_list = []
    for item in data:
        thing_idk(item, final_list, 0)

    with open(out_path, 'w') as file:
        json.dump(final_list, file, indent=2)


if __name__ == '__main__':
    start_time = time.time()
    main()
    end_time = time.time()
    log("Completed in", var=time_convert(end_time - start_time))
