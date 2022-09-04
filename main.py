from __future__ import annotations
import argparse
from src.train import train_model, evaluate_model, train_groups, eval_groups, extract_eval_data
from src.transforms import transform_coverage
from scripts.occurring_features import get_occurring_features, get_overlapping_features
from scripts.uniprotapi import multifasta_gen
import json

def main():
    choices = [
        'train',
        'evaluate',
        'extract',
        'features',
        'transform',
        'uniprot',
        'preprocess'
    ]
    parser = argparse.ArgumentParser()
    parser.add_argument('action', type=str, choices=['train', 'evaluate', 'extract', 'features', 'transform', 'uniprot'])
    parser.add_argument('--model', type=str)
    parser.add_argument('--gpath', type=str)
    parser.add_argument('--path', type=str, nargs='+')
    args = parser.parse_args()
    
    if args.action == 'train':
        if args.model == '-1':
            print("TEST")
            #train_groups()
            #eval_groups()
            extract_eval_data()
        elif args.model == 'all':
            train_groups()
        else:
            train_model(args.model)
    if args.action == 'evaluate':
        evaluate_model(args.model)
    if args.action == 'features':
        if len(args.path) == 1:
            print(get_occurring_features(args.path[0]))
        if len(args.path) == 2:
            print(get_overlapping_features(args.path[0], args.path[1]))
    if args.action == 'transform':
        annotation = json.load(open(args.path[0]))
        occurring_features = get_occurring_features(args.gpath)
        print(transform_coverage(annotation, occurring_features))
    if args.action == 'uniprot':
        multifasta_gen()
    if args.action == 'preprocess':
        pass

if __name__ == '__main__':
    main()
    
