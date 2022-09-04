import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json
from matplotlib.pyplot import figure
import matplotlib

def get_features(data, evalue_threshold=0, include_na_eval=False, cutoff=0.5):
    features = {}
    for protein in data.values():
        for feature, tools in protein.items():
            if tools and feature != 'length':
                for feature_type in tools:
                    feature_type = 'seg_LCR' if feature_type == 'seg_low_complexity_regions' else feature_type
                    if feature_type not in features:
                        features[feature_type] = 1
                    else:
                        features[feature_type] += 1
    #occurring_features = [feature for feature, count in features.items() if count >= len(data['feature'])//(1/cutoff)]
    #occurring_features = list(features.keys()) if not occurring_features else occurring_features
    features = {key: val/len(data) for key, val in features.items()}
    return features

def heatmap_gen(group):
    #figure(figsize=(8, 6), dpi=80)
    data = json.load(open(f'./data/anno/{group}.json'))['feature']
    print(len(data))
    prevalences = get_features(data)


    data_dict = {}
    proteins = []
    for p,anno in data.items():
        proteins.append(p)
        for tool, instances in anno.items():
            plength = anno['length']
            if tool != 'length':
                for domain, val in instances.items():
                    domain = 'seg_LCR' if domain == 'seg_low_complexity_regions' else domain
                    if domain not in data_dict:
                        data_dict[domain] = []
                    coverage = 0.0
                    prev_instance = None
                    for i in val['instance']:
                        if prev_instance:
                            if prev_instance[1] > i[0]:
                                coverage += (i[1] - prev_instance[1])/plength
                            else:
                                coverage += (i[1] - i[0])/plength
                        else:
                            coverage += (i[1] - i[0])/plength
                        prev_instance = i
                    data_dict[domain].append((p,coverage))
    
    for v in data_dict.values():
        add_to_v = []
        for p in proteins:
            if p not in [i[0] for i in v]:
                add_to_v.append((p,0.0))
        v += add_to_v
        v.sort(key=lambda x: x[0])
    
    
    df_dict = {k: [i[1] for i in v] for k,v in data_dict.items() if prevalences[k] >= 0.05}
    count = 0
    #for i in range(len(data_dict['pfam_IDH'])):
    #    if data_dict['pfam_IDH'][i][1] and data_dict['pfam_Iso_dh'][i][1]:
    #        count += 1
    #print("COUNT", count)
    df = pd.DataFrame(df_dict).T
    
    plt.imshow(df, cmap ="afmhot", interpolation='nearest', aspect='60', vmin=0.0, vmax=1.0)
    plt.colorbar(location='bottom', shrink=0.5)
    plt.xticks(fontsize=15)
    plt.yticks(range(len(df.index)), df.index, fontsize=13, rotation=60)
    plt.savefig(f'../{group}.png', dpi=1200)
    plt.show()

heatmap_gen('K00026')

def heatmap_gen_for_protein_set(pset, data_name):
    data = {pname: json.load(open(f'./data/eval/eval_data/{pname}.json')) for pname in pset}
    print(len(data))
    prevalences = get_features(data)


    data_dict = {}
    proteins = []
    for p,anno in data.items():
        proteins.append(p)
        for tool, instances in anno.items():
            plength = anno['length']
            if tool != 'length':
                for domain, val in instances.items():
                    domain = 'seg_LCR' if domain == 'seg_low_complexity_regions' else domain
                    if domain not in data_dict:
                        data_dict[domain] = []
                    coverage = 0.0
                    prev_instance = None
                    for i in val['instance']:
                        if prev_instance:
                            if prev_instance[1] > i[0]:
                                coverage += (i[1] - prev_instance[1])/plength
                            else:
                                coverage += (i[1] - i[0])/plength
                        else:
                            coverage += (i[1] - i[0])/plength
                        prev_instance = i
                    data_dict[domain].append((p,coverage))
    
    for v in data_dict.values():
        add_to_v = []
        for p in proteins:
            if p not in [i[0] for i in v]:
                add_to_v.append((p,0.0))
        v += add_to_v
        v.sort(key=lambda x: x[0])
    
    
    df_dict = {k: [i[1] for i in v] for k,v in data_dict.items() if prevalences[k] >= 0.01}
    df = pd.DataFrame(df_dict).T
    
    plt.imshow(df, cmap ="afmhot", interpolation='nearest', aspect='auto', vmin=0.0, vmax=1.0)
    plt.colorbar(location='bottom', shrink=0.5)
    plt.xticks(fontsize=15)
    plt.yticks(range(len(df.index)), df.index, fontsize=15, rotation=60)
    plt.savefig(f'../{data_name}.png', dpi=1200)
    plt.show()
    


k00241_fps = ['A0A0G2JNH3', 'A0A0U1RRN3', 'A0A0X1KG70', 'A0A126GWI2', 'A0A1B0GVL6', 'A0A1W2PQU2', 'A0A2R8Y5X9', 'A0A2R8YEG4', 
'A0A2R8YEH3', 'A0A2R8YEV3', 'A0A3B3IT45', 'A0PJK1', 'A0PK05', 'A1L3X0', 'A4D2G3', 'A6NC51', 'A6NCV1', 'A6NDH6', 
'A6NDL8', 'A6NDP7', 'A6NF89', 'A6NGC4', 'A6NH00', 'A6NHG9', 'A6NI61', 'A6NIJ9', 'A6NJY4', 'A6NJZ3', 'A6NKK0', 'A6NL08', 
'A6NL26', 'A6NL99', 'A6NM03', 'A6NM76', 'A6NML5', 'A6NMU1', 'A6NMZ5', 'A8MWL7', 'B0YJ81', 'B7ZAQ6', 'B9EJG8', 'G3V0H7', 
'H3BMG7', 'O00155', 'O00501', 'O00574', 'O14569', 'O14581', 'O14609', 'O14626', 'O14718', 'O14843', 'O14894', 'O15243', 
'O15374', 'O15375', 'O15403', 'O15427', 'O15529', 'O15551', 'O43194', 'O43246', 'O43603', 'O43688', 'O43749', 'O43826',
 'O43869', 'O60403', 'O60404', 'O60412', 'O60635', 'O60636', 'O60779', 'O75352', 'O76001', 'O76002', 'O76099', 'O76100',
  'O95006', 'O95007', 'O95013', 'O95047', 'O95136', 'O95214', 'O95221', 'O95222', 'O95371', 'O95377', 'O95452', 'O95471',
   'O95484', 'O95674', 'O95832', 'O95857', 'O95859', 'O95907', 'O95918', 'O95977', 'O95992', 'P00156', 'P00395', 'P00846',
    'P03891', 'P03901', 'P03915', 'P0C623', 'P0C628', 'P0C646', 'P0C7N5', 'P0C7T2', 'P0C7T3', 'P0CG08', 'P0DMS8', 'P0DMU2', 
    'P0DN80', 'P0DN81', 'P0DN82', 'P0DP42', 'P11168', 'P13945', 'P14416', 'P14672', 'P18577', 'P21462', 'P21731', 'P21926',
    'P25021', 'P25089', 'P25090', 'P25105', 'P28566', 'P29275', 'P29972', 'P30408', 'P30518', 'P30536', 'P30542', 'P30825', 
    'P30872', 'P31213', 'P31391', 'P32249', 'P32745', 'P34982', 'P35346', 'P35414', 'P41146', 'P41181', 'P41440', 'P43116', 
    'P43657', 'P46095', 'P46721', 'P47881', 'P47884', 'P48029', 'P48039', 'P51674', 'P52569', 'P53985', 'P54219', 'P56746', 
    'P56747', 'P56748', 'P56750', 'P56880', 'P57739', 'P58170', 'P58173', 'P58180', 'P58181', 'P58182', 'P59533', 'P59535', 
    'P59536', 'P59537', 'P59538', 'P59539', 'P59540', 'P59541', 'P59542', 'P59543', 'P59551', 'P59922', 'P60893', 'P78369', 
    'P78381', 'P78383', 'P82251', 'Q01362', 'Q02094', 'Q02161', 'Q03431', 'Q05940', 'Q0GE19', 'Q12908', 'Q13258', 'Q13286', 
    'Q13530', 'Q13571', 'Q13606', 'Q13607', 'Q13639', 'Q14330', 'Q14439', 'Q14656', 'Q14728', 'Q14973', 'Q15391', 'Q15612', 
    'Q15615', 'Q15617', 'Q15622', 'Q15629', 'Q16572', 'Q24JQ0', 'Q3KNW5', 'Q3LI62', 'Q495N2', 'Q49SQ1', 'Q53R12', 'Q53TN4', 
    'Q5GH77', 'Q5HYL7', 'Q5J8X5', 'Q5NUL3', 'Q5SR56', 'Q5TF39', 'Q5TGU0', 'Q5VZR4', 'Q5VZY2', 'Q5W0B7', 'Q63ZE4', 'Q695T7', 
    'Q6IEV9', 'Q6IEY1', 'Q6IF00', 'Q6IF42', 'Q6IF99', 'Q6IFG1', 'Q6IFH4', 'Q6N075', 'Q6NT16', 'Q6NTF9', 'Q6NUT3', 'Q6UW68', 
    'Q6UX65', 'Q6UXD7', 'Q6ZRR5', 'Q6ZVE7', 'Q71RS6', 'Q76EJ3', 'Q7RTP0', 'Q7RTR8', 'Q7RTY0', 'Q7RTY1', 'Q7Z2H8', 'Q7Z3Q1', 
    'Q7Z418', 'Q7Z5S9', 'Q7Z602', 'Q7Z7B1', 'Q86SM8', 'Q86SQ6', 'Q86TG1', 'Q86VZ1', 'Q86WB7', 'Q86WI0', 'Q8IXE1', 'Q8IZ96', 
    'Q8IZM9', 'Q8IZV2', 'Q8N0Y3', 'Q8N127', 'Q8N144', 'Q8N146', 'Q8N148', 'Q8N162', 'Q8N2M4', 'Q8N349', 'Q8N357', 'Q8N370', 
    'Q8N5M9', 'Q8N5U1', 'Q8N609', 'Q8N628', 'Q8N682', 'Q8N6F1', 'Q8N6L1', 'Q8N7C4', 'Q8N8Q1', 'Q8N8Q9', 'Q8NA29', 'Q8NBD8', 
    'Q8NBI2', 'Q8NBI5', 'Q8NBQ7', 'Q8NCC5', 'Q8NCK7', 'Q8NCS4', 'Q8NCU8', 'Q8NE00', 'Q8NFB2', 'Q8NFF2', 'Q8NG76', 'Q8NG77', 
    'Q8NG80', 'Q8NG85', 'Q8NG95', 'Q8NG97', 'Q8NG98', 'Q8NG99', 'Q8NGA0', 'Q8NGA1', 'Q8NGA2', 'Q8NGA5', 'Q8NGA6', 'Q8NGA8', 
    'Q8NGB4', 'Q8NGB6', 'Q8NGB8', 'Q8NGB9', 'Q8NGC2', 'Q8NGC3', 'Q8NGC5', 'Q8NGC8', 'Q8NGD0', 'Q8NGD1', 'Q8NGD2', 'Q8NGD3', 
    'Q8NGD5', 'Q8NGE0', 'Q8NGE1', 'Q8NGE2', 'Q8NGE3', 'Q8NGE5', 'Q8NGE8', 'Q8NGE9', 'Q8NGF3', 'Q8NGF6', 'Q8NGF7', 'Q8NGF8', 
    'Q8NGF9', 'Q8NGG0', 'Q8NGG1', 'Q8NGG4', 'Q8NGG8', 'Q8NGH8', 'Q8NGH9', 'Q8NGI2', 'Q8NGI4', 'Q8NGI6', 'Q8NGI7', 'Q8NGI8', 
    'Q8NGI9', 'Q8NGJ1', 'Q8NGJ3', 'Q8NGJ4', 'Q8NGJ5', 'Q8NGJ6', 'Q8NGK3', 'Q8NGK5', 'Q8NGL0', 'Q8NGL2', 'Q8NGL3', 'Q8NGL4', 
    'Q8NGL7', 'Q8NGL9', 'Q8NGM1', 'Q8NGM8', 'Q8NGM9', 'Q8NGN3', 'Q8NGN4', 'Q8NGN5', 'Q8NGN6', 'Q8NGP2', 'Q8NGP3', 'Q8NGP4', 
    'Q8NGP6', 'Q8NGP8', 'Q8NGP9', 'Q8NGQ2', 'Q8NGQ4', 'Q8NGQ5', 'Q8NGQ6', 'Q8NGR3', 'Q8NGS0', 'Q8NGS2', 'Q8NGS3', 'Q8NGS4', 
    'Q8NGS5', 'Q8NGS6', 'Q8NGS8', 'Q8NGS9', 'Q8NGT0', 'Q8NGT1', 'Q8NGT2', 'Q8NGT5', 'Q8NGT7', 'Q8NGT9', 'Q8NGU1', 'Q8NGU2', 
    'Q8NGV0', 'Q8NGX0', 'Q8NGX1', 'Q8NGX2', 'Q8NGX3', 'Q8NGX5', 'Q8NGX8', 'Q8NGX9', 'Q8NGY1', 'Q8NGY2', 'Q8NGY5', 'Q8NGY6', 
    'Q8NGY9', 'Q8NGZ0', 'Q8NGZ2', 'Q8NGZ3', 'Q8NGZ4', 'Q8NGZ6', 'Q8NH01', 'Q8NH04', 'Q8NH05', 'Q8NH08', 'Q8NH09', 'Q8NH16', 
    'Q8NH18', 'Q8NH19', 'Q8NH21', 'Q8NH37', 'Q8NH42', 'Q8NH43', 'Q8NH48', 'Q8NH50', 'Q8NH51', 'Q8NH53', 'Q8NH57', 'Q8NH60', 
    'Q8NH61', 'Q8NH64', 'Q8NH69', 'Q8NH70', 'Q8NH72', 'Q8NH73', 'Q8NH74', 'Q8NH79', 'Q8NH85', 'Q8NH90', 'Q8NH93', 'Q8NHA4', 
    'Q8NHA8', 'Q8NHC4', 'Q8NHC6', 'Q8NHC7', 'Q8NHC8', 'Q8NHS3', 'Q8TBG9', 'Q8TBR7', 'Q8TCB6', 'Q8TCT9', 'Q8TCU3', 'Q8TDB8', 
    'Q8TDU5', 'Q8TDU9', 'Q8TED4', 'Q8WUM9', 'Q8WUX1', 'Q8WW43', 'Q8WY07', 'Q8WZ84', 'Q8WZ92', 'Q8WZ94', 'Q8WZA6', 'Q92536', 
    'Q92633', 'Q92911', 'Q969K7', 'Q969L2', 'Q969S0', 'Q96B33', 'Q96B96', 'Q96BI1', 'Q96BI3', 'Q96CH1', 'Q96HV5', 'Q96IK0', 
    'Q96IZ2', 'Q96KK4', 'Q96LA9', 'Q96LB1', 'Q96LB2', 'Q96MC6', 'Q96MV1', 'Q96N19', 'Q96P67', 'Q96P88', 'Q96R08', 'Q96R09', 
    'Q96R28', 'Q96R30', 'Q96R45', 'Q96R47', 'Q96R48', 'Q96R54', 'Q96R69', 'Q96R72', 'Q96RA2', 'Q96RB7', 'Q96RC9', 'Q96RD0', 
    'Q96RD1', 'Q96RD2', 'Q96RD3', 'Q96RJ0', 'Q96S97', 'Q96SA4', 'Q96T55', 'Q99437', 'Q99463', 'Q99500', 'Q99679', 'Q9BRV3', 
    'Q9BRY0', 'Q9BVC6', 'Q9BVK8', 'Q9BW60', 'Q9BXA5', 'Q9BXC1', 'Q9BZD2', 'Q9BZW5', 'Q9GZK3', 'Q9GZK4', 'Q9GZK6', 'Q9GZM6', 
    'Q9GZP9', 'Q9H205', 'Q9H207', 'Q9H208', 'Q9H209', 'Q9H210', 'Q9H255', 'Q9H2C5', 'Q9H2V7', 'Q9H310', 'Q9H339', 'Q9H340', 
    'Q9H341', 'Q9H344', 'Q9H3U5', 'Q9H6F2', 'Q9H8P0', 'Q9HAB3', 'Q9HB14', 'Q9HBW0', 'Q9HC07', 'Q9HC97', 'Q9NP91', 'Q9NP94', 
    'Q9NPC1', 'Q9NPC2', 'Q9NPD5', 'Q9NPI0', 'Q9NQ40', 'Q9NQN1', 'Q9NQQ7', 'Q9NQS5', 'Q9NRA2', 'Q9NRX5', 'Q9NS66', 'Q9NS82', 
    'Q9NTQ9', 'Q9NUH8', 'Q9NUM3', 'Q9NVV5', 'Q9NWC5', 'Q9NWF4', 'Q9NY64', 'Q9NY91', 'Q9NYG8', 'Q9NYV7', 'Q9NYV8', 'Q9NYV9', 
    'Q9NYW0', 'Q9NYW1', 'Q9NYW2', 'Q9NYW3', 'Q9NYW4', 'Q9NYW5', 'Q9NYW6', 'Q9NYW7', 'Q9NZP0', 'Q9NZP2', 'Q9NZP5', 'Q9P055', 
    'Q9P0N5', 'Q9P0S9', 'Q9P1Q5', 'Q9UBR5', 'Q9UGF5', 'Q9UGF6', 'Q9UGF7', 'Q9UHI5', 'Q9UHI7', 'Q9UIG8', 'Q9UJ42', 'Q9UKL2', 
    'Q9UKP6', 'Q9UKR5', 'Q9UNL2', 'Q9UPQ8', 'Q9UPY5', 'Q9Y257', 'Q9Y289', 'Q9Y2T5', 'Q9Y2T6', 'Q9Y4A9', 'Q9Y585', 'Q9Y5P0', 
    'Q9Y672']

k00242_fps = ['A0A2R8Y5X9', 'P21145', 'Q3MUY2', 'A0A0B4J2F0', 'P58511', 'Q8NGA8', 'Q8NGV5', 'Q8N755', 'Q9GZW8', 'Q9H346', 
'P0C6T2', 'Q9Y342', 'Q9GZK4', 'Q86WK9', 'P51810', 'Q9NZD1', 'P41145', 'Q8NH76', 'Q96SL1', 'Q6ZSM3', 'Q8NBI2', 'Q12999', 
'P32248', 'Q6P499', 'P35410', 'P04201', 'P49683', 'O95562', 'P31358', 'V9GZ13', 'A0A1B0GUW7', 'P47881', 'Q8TDV2', 'A6NJY1',
 'Q8NGA0', 'O14494', 'Q16570', 'Q8NH21', 'Q8NH59', 'Q92535', 'Q969S6', 'Q8NGG3', 'O00631', 'Q9NZ42', 'Q9UHE8', 'Q9BY10', 
 'Q9GZN6', 'Q96R69', 'Q6ZVE7', 'Q8NGG2', 'Q8TBQ9', 'O15432', 'Q969W0', 'P0C646', 'Q8NGC1']


#heatmap_gen_for_protein_set(k00242_fps, 'k00242_fps')
