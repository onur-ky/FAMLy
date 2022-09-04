import os
import json


def get_features(data, evalue_threshold=0, include_na_eval=False, cutoff=0.5):
    features = {}
    for protein in data['feature'].values():
        for feature, tools in protein.items():
            if tools and feature != 'length':
                for feature_type in tools:
                    if feature_type not in features:
                        features[feature_type] = 1
                    else:
                        features[feature_type] += 1
    #occurring_features = [feature for feature, count in features.items() if count >= len(data['feature'])//(1/cutoff)]
    #occurring_features = list(features.keys()) if not occurring_features else occurring_features
    features = {key: val/len(data['feature']) for key, val in features.items()}
    return features


def main():
    anno_path = os.path.join(os.getcwd(), 'data', 'anno')

    for fname in os.listdir(anno_path):
        fpath = os.path.join(anno_path, fname)
        data = json.load(open(fpath, 'r'))
        features = get_features(data)

        outpath = os.path.join(os.getcwd(), 'data', 'metadata', f'{fname.split(".")[0]}.json')
        outfile = open(outpath, 'w')
        json.dump(features, outfile)


if __name__ == '__main__':
    main()
