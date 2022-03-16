'''
select_OG_all_sp.py
Script to capture Orthogroups shared by all species (NOT SC) from OrthoFinder output.
Created by Juan Manuel Lentijo Mondejar and Lisandra Benítez Álvarez on October 2021.
'''

import argparse
import csv
import textwrap

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Select Orthogroups with all species present (NOT SC)',
            formatter_class=argparse.RawDescriptionHelpFormatter,
      epilog=textwrap.dedent('''\
	for extract the orthogroups formed by all species and not SC  
      python select_OG_all_sp.py -f=Orthogroups.GeneCount.tsv -o=out.txt -c="species number"
      for extract the orthogroups formed by species (in 2,3,4 columns) and not SC  
      python select_OG_all_sp.py -f=Orthogroups.GeneCount.tsv -o=out.txt -c=2,3,4
      '''))
    parser.add_argument('-c', type=str, required=True, help='Columns to extract. The number of the columns of the samples that you want')
    parser.add_argument('-f', type=str, required=True, help='input file [./Orthogroups/Orthogroups.GeneCount.tsv file from Orthofinder output]')
    parser.add_argument('-o', type=str, required=True, help='output file')
    args = parser.parse_args()
    columns =args.c
    if "," in columns:
        columns = [int(c) for c in args.c.split(",")] 
    else:
        columns = [c for c in range(1,int(columns)+1)]
    output_samples = []

    with open(args.f, "r") as samples_file:
        tsvreader = list(csv.reader(samples_file, delimiter="\t"))
        output_samples.append([tsvreader[0][0]] + [tsvreader[0][col] for col in columns])
        for line in tsvreader[1:]:
            discard = False
            for pos, col in enumerate(line[1:], start=1):
                if int(col) <= 1 and pos in columns:
                    discard = True
                
            if discard is False:
                total = 0
                for value in [int(line[col]) for col in columns]:
                    total += value
                output_samples.append([line[0]] + [line[col] for col in columns] + [total])
    
    with open(args.o, "w+") as output_file:
        tsvwriter = csv.writer(output_file, delimiter="\t")
        tsvwriter.writerows(output_samples)

    print("{} selected genes.".format(len(output_samples) - 1))

