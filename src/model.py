import os
import torch
from torch import nn
from torch.utils.data import DataLoader, Subset, ConcatDataset
from collections import OrderedDict
from src.dataset import CustomAnnotationDataSet
from src.transforms import transform_coverage, transform_label
import logging
import numpy as np
import json
logging.getLogger().setLevel(logging.INFO)

DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'
LEARNING_RATE = 1e-3
BATCH_SIZE = 64
EPOCHS = 2
NUM_REGIONS = 30


GROUPS_TRAIN_PATH = os.path.join(
    os.getcwd(), 'data', 'train', 'groups_train'
)
GROUPS_TRAIN_LABELS = os.path.join(
    os.getcwd(), 'data', 'train', 'groups_train_labels.csv'
)


GROUPS_TEST_PATH = os.path.join(
    os.getcwd(), 'data', 'train', 'groups_train'
)
GROUPS_TEST_LABELS = os.path.join(
    os.getcwd(), 'data', 'train', 'groups_train_labels.csv'
)


COMP_TRAIN_PATH = os.path.join(
    os.getcwd(), 'data', 'train', 'comp_train'
)
COMP_TRAIN_LABELS = os.path.join(
    os.getcwd(), 'data', 'train', 'comp_train_labels.csv'
)


COMP_TEST_PATH = os.path.join(
    os.environ.get('ROOT_DIR'), 'data', 'train', 'comp_test'
)
COMP_TEST_LABELS = os.path.join(
    os.environ.get('ROOT_DIR'), 'data', 'train', 'comp_test_labels.csv'
)


def get_occurring_features(data):
    features = []
    for protein in data['feature']:
        for feature in data['feature'][protein]:
            if data['feature'][protein][feature] and feature != 'length':
                for feature_type in data['feature'][protein][feature]:
                    if feature_type not in features:
                        features.append(feature_type)
    return features


class LinearLayersModel(nn.Module):
    def __init__(self, input_nodes):

        super(LinearLayersModel, self).__init__()
        self.input_nodes = input_nodes

        self.flatten = nn.Flatten()
        self.linear_relu_stack = nn.Sequential(self.unpack_topology())
        self.linear_relu_stack.apply(self.init_weights)

    def forward(self, x):
        x = self.flatten(x)
        logits = self.linear_relu_stack(x)
        return logits

    def unpack_topology(self):
        """
        Takes the topology list and
        returns an ordered dict.
        """
        #self.topology[-1] += 1  # Two nodes in the last layer
        input_for_ordered_dict = []
        input_for_ordered_dict.append(
            (f'linear{1}', nn.Linear(self.input_nodes, self.input_nodes))
        )
        input_for_ordered_dict.append(
            (f'relu{1}', nn.ReLU())
        )
        input_for_ordered_dict.append(
            (f'linear{2}', nn.Linear(self.input_nodes, 1))
        )
        input_for_ordered_dict.append(
            (f'sigmoid{2}', nn.Sigmoid())
        )
        return OrderedDict(input_for_ordered_dict)
    
    @classmethod
    def init_weights(cls, m):
        if isinstance(m, nn.Linear):
            torch.nn.init.kaiming_uniform_(m.weight, mode='fan_in', nonlinearity='relu', a=0)
            m.bias.data.fill_(0.01)


class FAMLModel:
    def __init__(self, group_label):
        group_metadata_path = os.path.join(
            os.getcwd(), 'data', 'metadata', f'{group_label}.json'
        )
        eval_data_path = os.path.join(
            os.getcwd(), 'data', 'eval', 'eval_data'
        )
        self.eval_data_path = eval_data_path
        
        with open(group_metadata_path) as group_metadata_file:
            group_metadata = json.load(group_metadata_file)
            self.feature_space = [key for key, val in group_metadata.items() if val >= 0.50]
            input_nodes = sum([1 for val in group_metadata.values() if val >= 0.50])*NUM_REGIONS

        self.model = LinearLayersModel(input_nodes)
        self.group_label = group_label

        groups_train_dataset = CustomAnnotationDataSet(
            anno_dir=GROUPS_TRAIN_PATH,
            labels_file=GROUPS_TRAIN_LABELS,
            transform=transform_coverage,
            target_transform=transform_label,
            group=group_label
        )

        groups_test_dataset = CustomAnnotationDataSet(
            anno_dir=GROUPS_TEST_PATH,
            labels_file=GROUPS_TEST_LABELS,
            transform=transform_coverage,
            target_transform=transform_label,
            group=group_label
        )

        comp_train_dataset = CustomAnnotationDataSet(
            anno_dir=GROUPS_TEST_PATH,
            labels_file=GROUPS_TEST_LABELS,
            transform=transform_coverage,
            target_transform=transform_label,
            group=group_label
        )

        comp_test_dataset = CustomAnnotationDataSet(
            anno_dir=GROUPS_TEST_PATH,
            labels_file=GROUPS_TEST_LABELS,
            transform=transform_coverage,
            target_transform=transform_label,
            group=group_label
        )

        self.train_dataset = self.ds_gen(group_label, groups_train_dataset, comp_train_dataset)
        self.test_dataset = self.ds_gen(group_label, groups_test_dataset, comp_test_dataset)

        self.loss_fn = nn.BCEWithLogitsLoss()
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=1e-3)
        self.accuracy = None  # For tracking and visualization purposes
        self.loss = None
        self.weights_path = os.path.join(
            os.getcwd(), 'data', 'weights', f'{group_label}.pt'
        )

    @staticmethod
    def ds_gen(group_label, groups_ds, comp_ds):
        negative_comp_idx = np.where(comp_ds.targets != group_label)[0]
        negative_comp_subset = Subset(comp_ds, negative_comp_idx)

        negative_groups_idx = np.where(groups_ds.targets != group_label)[0]
        negative_groups_subset = Subset(groups_ds, negative_groups_idx)

        negative_ds = ConcatDataset([negative_comp_subset, negative_groups_subset])

        positive_groups_idx = np.where(groups_ds.targets == group_label)[0]
        positive_groups_subset = Subset(groups_ds, positive_groups_idx)

        positive_ds = ConcatDataset(
            [positive_groups_subset for i in
             range(len(negative_ds) // len(positive_groups_subset))]
        )

        ds = ConcatDataset([positive_ds, negative_ds])
        return ds

    def training_loop(self):
        dataloader = DataLoader(self.train_dataset, batch_size=BATCH_SIZE, shuffle=True)
        size = len(dataloader.dataset)

        for batch, (x, y, z) in enumerate(dataloader):
            prediction = self.model(x)
            loss = self.loss_fn(prediction, y)
            self.optimizer.zero_grad()
            loss.backward()
            self.optimizer.step()

            if batch % 100 == 0:
                loss, current = loss.item(), batch * len(x)
                logging.info(f'loss: {loss:>7f} [{current:>5d}/{size:>5d}]')

    def test_loop(self):
        dataloader = DataLoader(self.test_dataset, batch_size=BATCH_SIZE, shuffle=True)
        size = len(dataloader.dataset)
        num_batches = len(dataloader)
        test_loss, correct = 0, 0
        tps = 0
        tns = 0
        with torch.no_grad():
            for x, y, z in dataloader:
                prediction = self.model(x)
                #for i in range(len(prediction)):
                #    print(f'score: {prediction[i]} | {"positive" if y[i].item() == 1. else "negative"} | {z[i]}')
                #print("Above 0.9: ", len([i for i in prediction if i >= 0.9]))
                tp = len([prediction[i] for i in range(len(prediction)) if prediction[i] >= 0.9 and y[i].item() == 1.])
                tn = len([prediction[i] for i in range(len(prediction)) if prediction[i] <= 0.9 and y[i].item() == 0.])
                #print("True positives: ", tp)
                #print("True negatives: ", tn)
                #print("Accuracy: ", accuracy)
                #input("Enter to continue")
                test_loss += self.loss_fn(prediction, y).item()
                prediction = torch.tensor([[1.] if i[0] >= 0.9 else [0.] for i in prediction])
                test_result = torch.eq(prediction, y)
                correct += test_result.sum().item()
        test_loss /= num_batches
        correct /= size
        logging.info(f'Accuracy: {(100*correct):>0.1f}%, Avg loss: {test_loss:>8f}\n')
        return correct, test_loss

    def train(self):
        for i in range(EPOCHS):
            print(f'----------------------------\nEpoch: {i+1}')
            self.training_loop()
            accuracy, loss = self.test_loop()
            """
            if loss <= 0.6:
                self.accuracy = accuracy
                self.loss = loss
                return True
            """
        self.accuracy = accuracy
        self.loss = loss
        return False if accuracy < 0.95 else True # If false retrain
    
    def eval(self):
        outfile_path = os.path.join(
            os.getcwd(), 'data', 'eval', 'labels', f'{self.group_label}.csv'
        )
        self.load()
        with open(outfile_path, 'w+') as outfile:
            logging.info(f'Evaluating model {self.group_label}')
            for i, sample in enumerate(os.listdir(self.eval_data_path)):
                if not i%200:
                    logging.info(f'\t{i+1}/{len(os.listdir(self.eval_data_path))}')
                data = json.load(open(
                    os.path.join(self.eval_data_path, sample)
                ))
                transformed_data = transform_coverage(data, self.feature_space)
                prediction = self.model(transformed_data)
                #print(prediction)
                #input("press enter")
                outfile.write(f'{sample.split(".")[0]},{prediction[0][0].item()}\n')

    def save(self):
        model_data = {
            'model_state': self.model.state_dict(),
            'optimizer_state': self.optimizer.state_dict(),
            'loss': self.loss,
            'accuracy': self.accuracy
        }
        torch.save(model_data, self.weights_path)
        logging.info('Model data saved')
        #quit()
    
    def load(self):
        checkpoint = torch.load(self.weights_path)
        self.model.load_state_dict(checkpoint['model_state'])
        self.optimizer.load_state_dict(checkpoint['optimizer_state'])
        logging.info('Weights loaded')

"""
success = False
while not success:
    supermodel = FAMLModel(1)
    success = supermodel.train()
supermodel.save()
"""
