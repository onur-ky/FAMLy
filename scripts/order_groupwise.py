import argparse
import json
import os
from Bio import SeqIO
from utils import print_progress_bar


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', type=str)
    args = parser.parse_args()
    inpath = args.path
    print("Ordering by group...")
    for i, fname in enumerate(os.listdir(inpath)):
        print_progress_bar(i+1, len(os.listdir(inpath)))
        with open(os.path.join(inpath, fname)) as infile:
            data = SeqIO.parse(infile, "fasta")
            for record in data:
                gname = record.id.split('#')[0]
                gfile = open(os.path.join(os.getcwd(), 'data', 'fasta', f'{gname}.fasta'), 'a')
                gfile.write('>' + str(record.id) + '\n' + str(record.seq) + '\n')
                gfile.close()


if __name__ == '__main__':
    main()

