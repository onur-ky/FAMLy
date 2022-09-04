import json

def get_occurring_features_old(path):
    data = json.load(open(path))
    features = []
    for protein in data['feature']:
        for feature in data['feature'][protein]:
            if data['feature'][protein][feature] and feature != 'length':
                for feature_type in data['feature'][protein][feature]:
                    if feature_type not in features:
                        features.append(feature_type)
    return features

def get_occurring_features(path, evalue_threshold=0, include_NA_eval=False):
    data = json.load(open(path))
    features = {}
    for protein in data['feature'].values():
        for feature, tools in protein.items():
            if tools and feature != 'length':
                for feature_type in tools:
                    if feature_type not in features:
                        features[feature_type] = 1
                        #features.append(feature_type)
                    else:
                        features[feature_type] += 1

    occurring_features = [feature for feature, count in features.items() if count >= len(data['feature'])//2]
    for i,j in features.items():
        print(i, j)
    print("NUM SEQ:" ,len(data['feature']))
    occurring_features = list(features.keys()) if not occurring_features else occurring_features
    return occurring_features

def get_overlapping_features(path1, path2):
    feature_space_1 = get_occurring_features(path1)
    feature_space_2 = get_occurring_features(path2)
    return (set(feature_space_1) & set(feature_space_2))


