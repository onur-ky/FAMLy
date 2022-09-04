import os
import requests as r
import argparse
import json
from Bio import SeqIO


TRAIN_TO_TEST_CUTOFF = 0.1


def print_progress_bar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█', print_end='\r'):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end=print_end)
    # Print New Line on Complete
    if iteration == total:
        print()


def keggid_to_uprotid(keggid, species):
    url = f'https://rest.kegg.jp/conv/uniprot/{species}:{keggid}/'
    response = r.get(url)
    return response.text.split('\t')[1][3:-1]


def getid(header):
    return header.split('@')[1]


def species(entry):
    return entry.split('@')[0].split('#')[1]


def populate_comp(comp_proteome_path, comp_members_path):
    with open(comp_members_path) as comp_members_fasta, open(comp_proteome_path) as comp_proteome_anno:
        comp_members = SeqIO.parse(comp_members_fasta, 'fasta')
        comp_proteome_anno_data = json.load(comp_proteome_anno)['feature']

        print("\t\tGetting uniprot ids (this may take a while)...")
        comp_members = {keggid_to_uprotid(getid(i.id), species(i.id)): i.id.split('#')[0] for i in comp_members}

        comp_train_path = os.path.join(
            os.getcwd(), 'data', 'train', 'comp_train'
        )

        comp_test_path = os.path.join(
            os.getcwd(), 'data', 'train', 'comp_test'
        )

        comp_train_labels = open(os.path.join(
            os.getcwd(), 'data', 'train', 'comp_train_labels.csv'
        ), 'a+')

        comp_test_labels = open(os.path.join(
            os.getcwd(), 'data', 'train', 'comp_test_labels.csv'
        ), 'a+')

        print("\t\tPopulating complementary data directories")
        for i, (pid, annotation) in enumerate(comp_proteome_anno_data.items()):
            print_progress_bar(i+1, len(comp_proteome_anno_data))
            outdir = comp_test_path if i <= len(comp_proteome_anno_data)*TRAIN_TO_TEST_CUTOFF else comp_train_path
            labels = comp_test_labels if i <= len(comp_proteome_anno_data)*TRAIN_TO_TEST_CUTOFF else comp_train_labels
            gname = comp_members[pid] if pid in comp_members else 'other'

            labels.write(f'{pid},{gname}\n')

            with open(os.path.join(outdir, f"{pid}.json"), 'w+') as pannofile:
                json.dump(annotation, pannofile)

        comp_train_labels.close()
        comp_test_labels.close()


def populate_groups():
    cluster_path = os.path.join(
        os.getcwd(), 'data', 'cluster'
    )

    anno_path = os.path.join(
        os.getcwd(), 'data', 'anno'
    )

    groups_train_path = os.path.join(
        os.getcwd(), 'data', 'train', 'groups_train'
    )

    groups_test_path = os.path.join(
        os.getcwd(), 'data', 'train', 'groups_test'
    )

    groups_train_labels = open(os.path.join(os.getcwd(), 'data', 'train', 'groups_train_labels.csv'), 'a+')
    groups_test_labels = open(os.path.join(os.getcwd(), 'data', 'train', 'groups_test_labels.csv'), 'a+')

    repr_ids = []
    for i, fname in enumerate(os.listdir(cluster_path)):
        print_progress_bar(i+1, len(os.listdir(cluster_path)))
        with open(os.path.join(cluster_path, fname)) as gcluster:
            cluster_data = SeqIO.parse(gcluster, 'fasta')
            repr_ids += [p.id for p in cluster_data]

    print("\t\tPopulating complementary data directories")
    for fname in os.listdir(anno_path):
        anno_data = json.load(open(os.path.join(anno_path, fname)))['feature']
        print(f"\t\t\t{fname}")
        for i, (pid, annotation) in enumerate(anno_data.items()):
            print_progress_bar(i+1, len(anno_data))
            outdir = groups_test_path if i <= len(anno_data)*TRAIN_TO_TEST_CUTOFF else groups_train_path
            labels = groups_test_labels if i <= len(anno_data)*TRAIN_TO_TEST_CUTOFF else groups_train_labels
            gname = fname.split('.')[0]

            labels.write(f'{pid},{gname}\n')

            with open(os.path.join(outdir, f"{pid}.json"), 'w+') as pannodata:
                json.dump(annotation, pannodata)

    groups_train_labels.close()
    groups_test_labels.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--comp_gdata_path', type=str)
    parser.add_argument('--comp_whole_proteome_path', type=str)
    args = parser.parse_args()

    comp_proteome_path = args.comp_whole_proteome_path
    comp_members_path = args.comp_gdata_path
    print("\tPopulating complementary data...")
    populate_comp(comp_proteome_path, comp_members_path)
    print("\tPopulating group data...")
    populate_groups()


if __name__ == '__main__':
    main()
