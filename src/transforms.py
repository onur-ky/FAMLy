import torch
import os

def transform_coverage(annotation,
                       occurring_features,
                       evalue_threshold=999, 
                       include_NA_eval=True):
    num_of_regions = 30
    representation = {
        feature: {
            i+1: 0 for i in range(num_of_regions)
        } for feature in occurring_features
    }
    plength = annotation['length']
    region_length = max((plength // num_of_regions), 1)
    region_tailing_overhead_length = ((plength % num_of_regions) + region_length
                                        if plength > num_of_regions
                                        else 0)
    for feature in occurring_features:
        tool = 'coils2' if feature.split('_')[0] == 'coils' else feature.split('_')[0]
        if feature in annotation[tool]:

            evalue = float('NaN') if annotation[tool][feature]['evalue'] == 'NA' else annotation[tool][feature]['evalue']
            evalue_is_nan = torch.isnan(torch.tensor(evalue))

            if (not include_NA_eval and evalue_is_nan) or (evalue >= evalue_threshold):
                continue

            for instance in annotation[tool][feature]['instance']:
                instance_start_region = min(((instance[0]-1)//region_length)+1, num_of_regions)
                instance_end_region = min(((instance[1]-1)//region_length)+1, num_of_regions)

                if (instance_start_region == instance_end_region and 
                    instance_start_region == num_of_regions and 
                    region_tailing_overhead_length):
                    representation[feature][instance_start_region] = \
                        (instance[1]-instance[0]+1)/region_tailing_overhead_length

                elif (
                    (instance_start_region == instance_end_region and 
                    instance_start_region == num_of_regions and
                    not region_tailing_overhead_length) 
                    or
                    (instance_start_region == instance_end_region and
                    instance_start_region != num_of_regions)
                    ):
                    representation[feature][instance_start_region] = \
                        (instance[1]-instance[0]+1)/region_length
                
                else:
                    representation[feature][instance_start_region] = \
                        ((instance_start_region*region_length-instance[0])+1)/region_length
                    
                    for i in range(int(instance_start_region)+1, int(instance_end_region)):
                        representation[feature][i] = 1.0
                    
                    if instance_end_region == num_of_regions and region_tailing_overhead_length:
                        representation[feature][instance_end_region] = \
                            (instance[1]-((((instance_end_region-1)*region_length)+1))+1)/region_tailing_overhead_length
                    else:
                        representation[feature][instance_end_region] = \
                            (instance[1]-((((instance_end_region-1)*region_length)+1))+1)/region_length
    return torch.tensor([[[*pa.values()] for pa in representation.values()]], dtype=torch.float)


def transform_label(label, group_label):
    return torch.tensor([float(label == group_label)], dtype=torch.float)
