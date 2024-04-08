import os
import json
from tqdm import tqdm

dataset = 'Test'
adata_left_path = '/home/nuwaisir/Corridor/SPaSE/Data/King/Fixed_adatas/adata_Sham_1.h5ad'
adata_right_path = '/home/nuwaisir/Corridor/SPaSE/Data/King/Fixed_adatas/adata_Sham_1.h5ad'
adata_to_be_synthesized_path = adata_left_path
sample_left = 'Test_left'
sample_right = 'Test_right'
alpha = 0.001
lambda_sinkhorn = 0.01



sinkhorn = 1
dissimilarity = 'js'

mode = 1

numIterMaxEmd = 1000000
numInnerIterMax = 100
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

os.makedirs(os.path.dirname(config_path), exist_ok=True)

with open(f'{config_path}', 'w') as config_file:
    json.dump(config, config_file, indent=4)

os.system(f'python ../../../Workspace/SPaSE/main.py --config ../../../Workspace/SPaSE/configs/{config_file_name}')

with open(f'../../../Workspace/SPaSE/configs/{config_file_name}') as f:
    config = json.load(f)

config['mode'] = 2

with open(f'{config_path}', 'w') as config_file:
    json.dump(config, config_file, indent=4)

os.system(f'python ../../../Workspace/SPaSE/main.py --config ../../../Workspace/SPaSE/configs/{config_file_name}')
