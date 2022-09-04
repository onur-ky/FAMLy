import os
import torch
from torch.utils.data import Dataset, DataLoader
import pandas as pd
import json
from dotenv import load_dotenv
from src.transforms import transform_coverage, transform_label

load_dotenv()


class CustomAnnotationDataSet(Dataset):
    def __init__(self, labels_file, anno_dir, 
                 transform, target_transform, group):
        self.anno_labels = pd.read_csv(labels_file)
        self.anno_dir = anno_dir
        self.transform = transform
        self.target_transform = target_transform
        self.group = group

    def __len__(self):
        return len(self.anno_labels)
    
    def __getitem__(self, idx):
        anno_name = self.anno_labels.iloc[idx, 0]
        anno_path = os.path.join(self.anno_dir, f'{anno_name}.json')
        label = self.anno_labels.iloc[idx, 1]

        group_metadata_path = os.path.join(
            os.getcwd(), 'data', 'metadata', f'{self.group}.json'
        )
        occurring_features = json.load(open(group_metadata_path))
        occurring_features = [domain for domain, prevalence in occurring_features.items() if prevalence >= 0.50]
        #print("OCCURRING FEATURES: ", occurring_features)
        with open(anno_path, 'r') as anno_file:
            annotation = json.load(anno_file)
        if self.transform:
            annotation = self.transform(annotation, occurring_features, include_NA_eval=True)
        if self.target_transform:
            label = self.target_transform(label, self.group)
        #print(annotation, anno_name)
        #print("-----------------------")
        return annotation, label, anno_name

    @property
    def targets(self):
        return self.anno_labels.iloc[:, 1].to_numpy()
