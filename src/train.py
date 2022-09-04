import os
import json
import logging

from matplotlib.pyplot import table
from src.model import FAMLModel
from dotenv import load_dotenv
logging.getLogger().setLevel(logging.INFO)
load_dotenv()

def train_groups():
    exclude = ['K00024', 'K00030', 'K00031', 'K00627']
    metadata_path = os.path.join(
        os.getcwd(), 'data', 'metadata'
    )
    for fname in os.listdir(metadata_path):
        gname = fname.split('.')[0]
        if gname not in exclude:
            logging.info(f'Training a model for the group with label {gname}')
            model = FAMLModel(gname)
            model.train()
            model.save()


def train_model(model_label):
    logging.info(f'Training a model for the group with label {model_label}')
    end_training = False
    while not end_training:
        model = FAMLModel(model_label)
        model.train()
        end_training = input("Retrain? [y/n]: ") == "n"
    model.save()


def evaluate_model(model_label):
    model = FAMLModel(model_label)
    model.eval()
    human_tca_list = {
        'K01647': ['O75390'],
        'K01648': ['P53396'],
        'K01681': ['Q99798', 'P21399'],
        'K00031': ['O75874', 'P48735'],
        'K00030': ['O43837', 'P51553', 'P50213'],
        'K00164': ['Q9ULD0', 'Q02218'],
        'K00658': ['P36957'],
        'K00382': ['P09622'],
        'K01899': ['P53597'],
        'K01900': ['Q96I99', 'Q9P2R7'],
        'K00234': ['P31040'],
        'K00235': ['P21912'],
        'K00236': ['Q99643'],
        'K00237': ['O14521'],
        'K01679': ['P07954'],
        'K00025': ['P40925'],
        'K00026': ['P40926'],
        'K01958': ['P11498'],
        'K01596': ['Q16822', 'P35558'],
        'K00161': ['P29803', 'P08559'],
        'K00162': ['P11177'],
        'K00627': ['P10515'],
        'K13997': ['O00330'],
    }
    eval_labels_path = os.path.join(
        os.getcwd(), 'data', 'eval', 'labels', f'{model_label}.csv'
    )
    eval_labels = open(eval_labels_path, 'r').read().split('\n')
    positives = human_tca_list[model_label] if model_label in human_tca_list else []
    true_positives = []
    false_positives = []
    results_dict = {
        'true_positives': [],
        'false_positives': [],
        'false_negatives': []
    }
    for sample in eval_labels[:-1]:
        if float(sample.split(',')[1]) > 0.99:
            if sample.split(',')[0] in positives:
                true_positives.append(sample.split(',')[0])
                results_dict['true_positives'].append(sample.split(',')[0])
            else:
                false_positives.append(sample.split(',')[0])
                results_dict['false_positives'].append(sample.split(',')[0])
        elif sample.split(',')[0] in positives:
            results_dict['false_negatives'].append(sample.split(',')[0])
    print('TRUE POSITIVES: \n', len(true_positives))
    print('FALSE POSITIVES: \n', len(false_positives))
    print('FALSE NEGATIVES: \n', len(positives)-len(true_positives), results_dict['false_negatives'])


def eval_groups():
    metadata_path = os.path.join(
        os.getcwd(), 'data', 'metadata'
    )
    for fname in os.listdir(metadata_path):
        gname = fname.split('.')[0]
        model = FAMLModel(gname)
        model.eval()

def extract_eval_data():
    human_tca_list = {
        'K01647': ['O75390'],
        'K01648': ['P53396'],
        'K01681': ['Q99798', 'P21399'],
        'K00031': ['O75874', 'P48735'],
        'K00030': ['O43837', 'P51553', 'P50213'],
        'K00164': ['Q9ULD0', 'Q02218'],
        'K00658': ['P36957'],
        'K00382': ['P09622'],
        'K01899': ['P53597'],
        'K01900': ['Q96I99', 'Q9P2R7'],
        'K00234': ['P31040'],
        'K00235': ['P21912'],
        'K00236': ['Q99643'],
        'K00237': ['O14521'],
        'K01679': ['P07954'],
        'K00025': ['P40925'],
        'K00026': ['P40926'],
        'K01958': ['P11498'],
        'K01596': ['Q16822', 'P35558'],
        'K00161': ['P29803', 'P08559'],
        'K00162': ['P11177'],
        'K00627': ['P10515'],
        'K13997': ['O00330'],
    }
    no_human_member = ['K00174', 'K00116', 'K00171', 'K03737', 'K00024', 'K00169', 'K18859', 'K00240', 'K00242', 
    'K05942', 'K00239', 'K00170', 'K01610', 'K15231', 'K01678', 'K15230', 'K18860', 'K17753', 'K01960', 'K01682',
     'K01959', 'K00177', 'K01676', 'K00246', 'K18118', 'K00163', 'K00172', 'K00241', 'K00176', 'K00175', 'K01903', 
     'K00247', 'K01902', 'K00245', 'K01677', 'K00244']
    metadata_path = os.path.join(
        os.getcwd(), 'data', 'metadata'
    )

    hits = []
    fps = []
    all_accuracy = []
    all_precision = []

    tablestrings = []
    for fname in os.listdir(metadata_path):
        gname = fname.split('.')[0]
        #print(f'--------------GROUP: {gname}--------------')
        eval_labels_path = os.path.join(
            os.getcwd(), 'data', 'eval', 'labels', f'{gname}.csv'
        )
        eval_labels = open(eval_labels_path, 'r').read().split('\n')
        positives = human_tca_list[gname] if gname in human_tca_list else []
        true_positives = []
        false_positives = []
        results_dict = {
            'true_positives': [],
            'false_positives': [],
            'false_negatives': []
        }
        for sample in eval_labels[:-1]:

            if float(sample.split(',')[1]) > 0.99:
                if sample.split(',')[0] in positives:
                    true_positives.append(sample.split(',')[0])
                    results_dict['true_positives'].append(sample.split(',')[0])
                else:
                    false_positives.append(sample.split(',')[0])
                    results_dict['false_positives'].append(sample.split(',')[0])
            elif sample.split(',')[0] in positives:
                results_dict['false_negatives'].append(sample.split(',')[0])
        #print('TRUE POSITIVES: \n', len(true_positives))
        #print('FALSE POSITIVES: \n', len(false_positives))
        #print('FALSE NEGATIVES: \n', len(positives)-len(true_positives))
        false_negatives = len(positives)-len(true_positives)
        true_negatives = 20874 - len(true_positives) + len(false_positives) + (len(positives)-len(true_positives))
        #print('TRUE NEGATIVES: \n', 20874 - len(true_positives) + len(false_positives) + (len(positives)-len(true_positives)))

        hits += true_positives
        fps += false_positives
        #print(f'{group}: {len(true_positives)},{len(false_positives)},{len(positives)-len(true_positives)}')
        true_positives = len(true_positives)
        false_positives = len(false_positives)
        accuracy = (true_positives + true_negatives) / (true_positives + true_negatives + false_positives + false_negatives)
        #print("ACCURACY: ", accuracy)
        precision = true_positives / (true_positives + false_positives) if true_positives or false_positives else 1.0
        #print("PRECISION: ", precision)
        #print("RECALL: ", true_positives / len(positives) if positives else 1.0)
        tableformatstring = f'{gname} & {true_positives} & {false_positives} & {false_negatives} & {("%.2f" % ((true_positives / len(positives))*100)) + "%" if positives else "-"} & {("%.2f" % (accuracy*100)) + "%"} & {("%.2f" % (precision*100)) + "%"} \\\\ \n\\hline'
        tableformatstring = tableformatstring.replace('%', '\\%')
        if gname not in no_human_member and (true_positives or false_positives or false_negatives):
            tablestrings.append(tableformatstring)
        all_accuracy.append(accuracy)
        all_precision.append(precision)
        outpath = os.path.join(
            os.getcwd(), 'data', 'eval', 'results', f'{gname}.json'
        )
        outfile = open(outpath, 'w+')
        json.dump(results_dict, outfile)
    print(hits, len(fps))
    print("AVERAGE ACCURACY: ", sum(all_accuracy)/len(all_accuracy))
    print("AVERAGE PRECISION: ", sum(all_precision)/len(all_precision))
    tablestrings = sorted(tablestrings)
    for _ in tablestrings:
        print(_)


#train_groups()
#eval_groups()
#extract_eval_data()

#model = FAMLModel("K00031")
