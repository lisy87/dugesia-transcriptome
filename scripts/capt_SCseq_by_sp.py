'''
capt_SCseq_by_sp.py
Script to capture the the name of Single Copy sequences by sample
Created by Juan Manuel Lentijo Mondejar and Lisandra Benítez Álvarez on June 2021.
'''

import argparse
import csv
from collections import defaultdict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Capture the name of Single Copy sequences by sample')
    parser.add_argument('-f', type=str, help='list of species [species_list.txt]')
    parser.add_argument('-i', type=str, help=' Single Copy identifiers by species [SingCopy_ID.txt]')
    args = parser.parse_args()
    data = defaultdict(list)

    with open(args.f, "r") as species_list:
        filenames = species_list.readlines()
        species_list.close()

    with open(args.i, "r") as sc_ids:
        csvreader = csv.reader(sc_ids, delimiter="\t")
        for row in csvreader:
            for (i, v) in enumerate(row):
                if i > 0:
                    data[i-1].append(v)
        sc_ids.close()

    for (i, v) in data.items():
        with open(filenames[i].strip(), "w+") as output:
            output.writelines("%s\n" % l for l in v)
            output.close()
            
