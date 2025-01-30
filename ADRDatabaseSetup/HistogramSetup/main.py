import json
import time
from decimal import Decimal

"""
{
name: Allele
type: 'bar'
stack: 'total'
data: [Freq]
}
"""

filename = "histogram.json"

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


def main():
    final_list = []
    for item in data:
        if item["Freq"] > 0:
            drug_stuff = item["Drugs"].strip().split(" | ")
            final_list.append({
                "name": item["Allele"],
                "type": "bar",
                "stack": 'bar',
                "data": [move_decimal_right(item["Freq"])],
                "loci": item["Gene"],
                "region": item["Pop"],
                "drug": drug_stuff})

    with open(filename + ".out.json", 'w') as file:
        json.dump(final_list, file, indent=2)


if __name__ == '__main__':
    begin_time = time.time()
    main()
    print("Total time taken : ", time.time()-begin_time)
