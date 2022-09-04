import os
from dotenv import load_dotenv
import logging
import json
import random  # For random sampling of training/test datasets

load_dotenv()
logging.getLogger().setLevel(logging.INFO)

class Preprocesser:
    def __init__(self):
        data = self.fetch_data()
        self.labels = data['labels']
        
    
    @classmethod
    def fetch_data(cls):
        groups_to_labels_file = os.path.join(
            os.environ.get('ROOT_DIR'), 'data', 'groups_to_labels.json'
        )
        labels_dict = json.load(open(groups_to_labels_file))
        data = {
            'labels': labels_dict
        }
        return data


    def seperate_raw_data(self, raw_data_path, train_data_path, test_data_path):
        for i, fname in enumerate(os.listdir(raw_data_path)):

            logging.info(f'Processing {fname}:\
                {i+1}/{len(os.listdir(raw_data_path))}')

            suffix = '' if fname != 'YEAST@559292@3.json' else 'YEAST#'

            group_file = open(os.path.join(raw_data_path, fname))
            group_data = json.load(group_file)
            annotations = group_data['feature']

            switch_after = len(annotations.items()) // 10
            switch_after = 1 if switch_after == 0 else switch_after
            for i, (protein_name, annotation) in enumerate(annotations.items()):
                if i < switch_after:
                    outfile_path = os.path.join(test_data_path, f'{suffix}{protein_name}.json')
                else:
                    outfile_path = os.path.join(train_data_path, f'{suffix}{protein_name}.json')
                with open(outfile_path, 'w+') as outfile:
                    json.dump(annotation, outfile)


    def populate_test_data(self, train_data_path, test_data_path):

        data_files = os.listdir(train_data_path)
        num_of_files_to_move = len(data_files) // 10  # 0.1 ratio
        
        for i in range(num_of_files_to_move):
            if not i % 100: 
                logging.info(f'Moving 0.1 of the files to test:\
    {i}/{num_of_files_to_move} done.')
            chosen_file = random.choice(data_files)
            os.rename(
                os.path.join(train_data_path, chosen_file),
                os.path.join(test_data_path, chosen_file)
            )
            data_files = os.listdir(train_data_path)  # Update list of files so we don't get duplicates


    def generate_labels(self, path_to_dataset, outfile_name):
        outfile_path = os.path.join(
            os.environ.get('ROOT_DIR'), 'data', 'train', outfile_name
        )
        with open(outfile_path, 'w+') as outfile:
            for i, fname in enumerate(os.listdir(path_to_dataset)):
                if not i % 100:
                    logging.info(f'Generating labels for dataset:\
    \n{path_to_dataset}.\
    \n{i}/{len(os.listdir(path_to_dataset))} done.')
                outfile.write(f'{fname},{self.labels[fname.split("#")[0]]}\n')


    def preprocess_raw_data(self):
        data_path = os.path.join(os.environ.get('ROOT_DIR'), 'data')
        raw_data_path = os.path.join(data_path, 'train', 'raw')
        train_data_path = os.path.join(data_path, 'train', 'train_data')
        test_data_path = os.path.join(data_path, 'train', 'test_data')

        self.seperate_raw_data(raw_data_path, train_data_path, test_data_path)
        #self.populate_test_data(train_data_path, test_data_path)
        self.generate_labels(train_data_path, 'train_labels.csv')
        self.generate_labels(test_data_path, 'test_labels.csv')

    @classmethod
    def preprocess_raw_eval_data(cls, fname):
        eval_raw_data_path = os.path.join(
            os.environ.get('ROOT_DIR'), 'data', 'eval', 'raw', fname
        )
        eval_data_path = os.path.join(
            os.environ.get('ROOT_DIR'), 'data', 'eval', 'eval_data'
        )

        logging.info(f'Processing evaluation data {fname}')
        data = json.load(open(eval_raw_data_path))
        for protein, annotation in data['feature'].items():
            outfile_path = os.path.join(eval_data_path, f'{protein}.json')
            with open(outfile_path, 'w+') as outfile:
                json.dump(annotation, outfile)

p = Preprocesser()
#p.preprocess_raw_eval_data('HUMAN@9606@3.json')