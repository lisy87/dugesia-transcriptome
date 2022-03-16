'''
sel_sc_by_sample.py
Script to capture Single Copy Orthogroups shared only by a group of species from OrthoFinder output.
Created by Juan Manuel Lentijo Mondejar and Lisandra Benítez Álvarez on June 2021.
'''
import argparse
import csv
import textwrap

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Select single copy genes by species using the Orthogroups.tsv file from Orthofinder output',
	    formatter_class=argparse.RawDescriptionHelpFormatter,
      epilog=textwrap.dedent('''\
      for extract the sc for the species in columns 1, 2, and 3 
      python sel_sc_by_sample.py -f=Orthogroups.tsv -o=out.txt -c=1,2,3
      '''))
    parser.add_argument('-c', type=str, required=True, help='columns to extract. The number of the columns of the samples that you want to extract starting by 1')
    parser.add_argument('-f', type=str, required=True, help='input file [Orthogroups.tsv from orthofinder output]')
    parser.add_argument('-o', type=str, required=True, help='output file')
    args = parser.parse_args()
    columns = [int(c) for c in args.c.split(",")]
    output_samples = []

    with open(args.f, "r") as samples_file:
        tsvreader = list(csv.reader(samples_file, delimiter="\t"))
        output_samples.append([tsvreader[0][0]] + [tsvreader[0][col] for col in columns])
        for line in tsvreader[1:]:
            discard = False
            for pos, col in enumerate(line[1:], start=1):
                if col != "" and not pos in columns:
                    discard = True
                if col == "" and pos in columns:
                    discard = True
                if "," in col and pos in columns:
                    discard = True

            if discard is False:
                output_samples.append([line[0]] + [line[col] for col in columns])
    
    with open(args.o, "w+") as output_file:
        tsvwriter = csv.writer(output_file, delimiter="\t")
        tsvwriter.writerows(output_samples)

    print("{} selected genes.".format(len(output_samples) - 1))
    
