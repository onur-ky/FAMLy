import requests
import os
import json

def get_fps():
    fps = set()
    eval_results_path = os.path.join(
        os.getcwd(), 'data', 'eval', 'eval_results', '99'
    )
    for fname in os.listdir(eval_results_path):
        data = json.load(
            open(os.path.join(eval_results_path, fname))
        )
        for fp in data['false_positives']:
            fps.add(fp)
    return fps

def get_fasta(uniprotid):
    url = f"http://www.uniprot.org/uniprot/{uniprotid}.fasta"
    response = requests.post(url)
    fasta = response.text
    return fasta

def multifasta_gen():
    multifasta = ''
    uniprotids = get_fps()
    for id in uniprotids:
        fasta = get_fasta(id)
        multifasta += fasta
    with open('fps.fasta', 'w+') as outfile:
        outfile.write(multifasta)
    print("Fasta generated")
