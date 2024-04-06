import os
import json
from tqdm import tqdm

dataset = 'None'
adata_left_path = 'None'
adata_right_path = 'None'
adata_to_be_synthesized_path = 'None'

cases = []

datasets = ['King_fixed']

sample_pairs = {
    'King_fixed': [['Sham_1', 'D3_1'], ['Sham_1', 'D3_3'], ['Sham_1', 'D7_2'], ['Sham_1', 'D7_3']],
}

sample_alpha_map = {
    '1hr': 0.01,
    '4hr': 0.01,
    'D3_1': 0.0001,
    'D3_3': 0.01,
    'D7_2': 0.0001,
    'D7_3': 0.01,
}
sample_lambda_map = {
    '1hr': 0.1,
    '4hr': 0.01,
    'D3_1': 0.1,
    'D3_3': 0.1,
    'D7_2': 0.1,
    'D7_3': 0.01,
}

dissimilarities = ['js']
sinkhorn_options = [1]

for dataset in datasets:
    for sample_pair in sample_pairs[dataset]:
        sample_left = sample_pair[0]
        sample_right = sample_pair[1]
        for dissimilarity in dissimilarities:
            alpha_options = [sample_alpha_map[sample_right]]
            lambda_options = [sample_lambda_map[sample_right]]
            for sinkhorn in sinkhorn_options:
                if sinkhorn == 1:
                    for alpha in alpha_options:
                        for lambda_sinkhorn in lambda_options:
                            cases.append({
                                'dataset': dataset,
                                'sample_left': sample_left,
                                'sample_right': sample_right,
                                'dissimilarity': dissimilarity,
                                'sinkhorn': sinkhorn,
                                'alpha': alpha,
                                'lambda_sinkhorn': lambda_sinkhorn
                            })
                else:
                    for alpha in alpha_options:
                        cases.append({
                            'dataset': dataset,
                            'sample_left': sample_left,
                            'sample_right': sample_right,
                            'dissimilarity': dissimilarity,
                            'sinkhorn': sinkhorn,
                            'alpha': alpha,
                            'lambda_sinkhorn': 1
                        })

for case in tqdm(cases):
    mode = 1
    dataset = case['dataset']
    sample_left = case['sample_left']
    sample_right = case['sample_right']
    lambda_sinkhorn = case['lambda_sinkhorn']
    sinkhorn = case['sinkhorn']
    dissimilarity = case['dissimilarity']
    alpha = case['alpha']
    numIterMaxEmd = 1000000
    numInnerIterMax = 10000
    init_map_scheme = "uniform"
    use_gpu = 1
    QC = 0

    config = {
        "mode": mode,
        "dataset": dataset,
        "sample_left": sample_left,
        "sample_right": sample_right,
        "adata_left_path": adata_left_path,
        "adata_right_path": adata_right_path,
        "adata_to_be_synthesized_path": adata_to_be_synthesized_path,
        "sinkhorn": sinkhorn,
        "lambda_sinkhorn": lambda_sinkhorn,
        "dissimilarity": dissimilarity,
        "alpha": alpha,
        "init_map_scheme": init_map_scheme,
        "numIterMaxEmd": numIterMaxEmd,
        "numInnerIterMax": numInnerIterMax,
        "use_gpu": use_gpu,
        "QC": QC,
        "data_folder_path": "../../../Data",
        "sample_left_hvg_h5_save_path": "../../../Data/King/Preprocessed",
        "sample_right_hvg_h5_save_path": "../../../Data/King/Preprocessed",
        "results_path": "../../../Workspace/SPaSE/results",
        "grid_search": 0,
    }

    config_file_name = f'config_{dataset}_{sample_left}_vs_{sample_right}_{dissimilarity}'
    if sinkhorn: config_file_name += f'_sinkhorn_lambda_{lambda_sinkhorn}_alpha_{alpha}.json'
    else:
        config_file_name += f'_alpha_{alpha}.json'

    config_path = f'../../../Workspace/SPaSE/configs/{config_file_name}'
    with open(f'{config_path}', 'w') as config_file:
        json.dump(config, config_file, indent=4)

    os.system(f'python ../../../Workspace/SPaSE/main.py --config ../../../Workspace/SPaSE/configs/{config_file_name}')

    with open(f'../../../Workspace/SPaSE/configs/{config_file_name}') as f:
        config = json.load(f)

    config['mode'] = 2

    with open(f'{config_path}', 'w') as config_file:
        json.dump(config, config_file, indent=4)

    os.system(f'python ../../../Workspace/SPaSE/main.py --config ../../../Workspace/SPaSE/configs/{config_file_name}')
