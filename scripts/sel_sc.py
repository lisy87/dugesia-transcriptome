'''
 sel_sc.py
Script to capture the nucleotide sequence of Single Copy genes from OrthoFinder and Transdecoder output.
Created by Juan Manuel Lentijo Mondejar and Lisandra Benítez Álvarez on November 2020.
'''

import gzip
import os
import argparse
import glob
import time
import textwrap


def read_file(filename):
    with open(filename, "r") as f:
        lines = f.read().splitlines()
        f.close()
    return lines


def read_database(filename):
    with open(filename, "r") as f:
        database = {}
        line = f.readline()
        while line:
            database[line.strip()] = f.readline().strip()
            line = f.readline()
        f.close()
    return database


def select_single_copy(filename, single_copy, database):
    select_single_copy = []
    with open(filename, "w") as f:
        for name in single_copy:
            name_gen = ">{}".format(name.strip())
            f.write(name_gen + "\n")
            f.write(database[name_gen] + "\n")
            select_single_copy.append(name_gen)
            select_single_copy.append(database[name_gen])
    f.close()
    return select_single_copy


def create_all_files(root_dir, all_seq_sel):
    filenames = os.path.join(root_dir, "Orthogroups_SingleCopyOrthologues.txt")
    with open(filenames, "r") as f:
        output_file = f.read().splitlines()
    f.close()
    line = 0
    for name in output_file:
        with open(os.path.join(root_dir, name + ".fasta"), "w") as f:
            for species_name, gen in all_seq_sel.items():
                f.write(">" + species_name + "_" + gen[line][1:] + "\n")
                f.write(gen[line + 1] + "\n")
        line += 2


def main(root_dir="./"):
    txt_files = [f.split(os.path.sep)[-1] for f in glob.glob(os.path.join(root_dir, "*seq.txt"))]
    all_seq_sel = {}
    for txt_name in txt_files:
        species_name = txt_name.replace("_seq.txt","")
        fasta_name = species_name + ".fasta"
        single_copy = read_file(filename=os.path.join(root_dir, txt_name))
        database = read_database(filename=os.path.join(root_dir, fasta_name))
        filename = os.path.join(root_dir, fasta_name.split(".")[0] + "_sc.fasta")
        all_seq_sel[species_name] = select_single_copy(filename, single_copy, database)
    create_all_files(root_dir, all_seq_sel)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Select single copy',
        formatter_class=argparse.RawDescriptionHelpFormatter,
      epilog=textwrap.dedent('''\
        IMPORTANT INFORMATION FOR THE USER!!!
        Put into a folder:
        1. The sequential CDS files where you go to look for the SC genes sequence "sample.fasta"
        2. The "sample_seq.txt" files, where you have the SC sequence name for every sample 
        3. The orthofinder file output "Orthogroups_SingleCopyOrthologues.txt"
        IMPORTANT NOTE!!!! You have to keep this structure in the names of the files
        Then, you can run the script IN THIS DIRECTORY just as: python sel_sc_edit.py
        If you use the -d option you can define an output  directory otherwise, the results will be write in the default directory "./"
        After run the script you get:
          - a fasta file for each of the species with all the OGs "sample_sc.fasta"
          - a fasta file for each of the OGs with all the species, ready to be aligned "OG*.fasta"
        '''))
    parser.add_argument('-d', type=str, help='output directory', default="." + os.path.sep)
    args = parser.parse_args()
    total_time = time.time()
    main(root_dir=args.d)
    print("Time elapsed in: {} seconds.".format(time.time() - total_time))
