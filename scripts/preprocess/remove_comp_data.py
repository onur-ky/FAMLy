import json
import os
from Bio import SeqIO
import argparse


def print_progress_bar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█', print_end='\r'):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end=print_end)
    # Print New Line on Complete
    if iteration == total:
        print()


def species(entry):
    return entry.split('@')[0].split('#')[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--remove', type=str)
    args = parser.parse_args()

    annopath = os.path.join(
        os.getcwd(), 'data', 'anno'
    )

    remove = args.remove

    for i, fname in enumerate(os.listdir(annopath)):
        print_progress_bar(i+1, len(os.listdir(annopath)))
        with open(os.path.join(annopath, fname), 'r+') as infile:
            data = json.load(infile)['feature']
            data = {key: val for key, val in data.items() if species(key) != remove}
            infile.seek(0)
            infile.truncate()
            json.dump(data, infile)


if __name__ == '__main__':
    main()
