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

def main():
    human_proteome_path = '/home/onur/Desktop/HUMAN@9606@3.json'
    human_data = json.load(open(human_proteome_path))['feature']
    eval_data_path = '/home/onur/src/FAML/data/eval/eval_data'
    for i, (key, val) in enumerate(human_data.items()):
        print_progress_bar(i+1, len(human_data))
        with open(f'{eval_data_path}/{key}.json', 'w+') as outfile:
            json.dump(val, outfile)

main()