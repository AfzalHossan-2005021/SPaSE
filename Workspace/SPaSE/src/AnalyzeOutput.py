import numpy as np
import pandas as pd
import torch
import ot
import scipy
import matplotlib.pyplot as plt
import warnings
from .DataLoader import DataLoader
from .Preprocessor import Preprocessor
import os
from .utils import plot_slice_pairwise_alignment_modified, calculate_cost_matrix, mirror, rotate, get_2hop_adatas, compute_null_distribution, visualize_goodness_of_mapping, scale_coords, QC, paste_pairwise_align_modified
from scipy.stats import variation
import scanpy as sc
import json
from scipy import stats
from sklearn import metrics
import seaborn as sns
from scipy.stats import gamma

class AnalyzeOutput:
    def __init__(self, config):
        self.config = config
        self.dataset = config['dataset']
        self.sample_left = config['sample_left']
        self.sample_right = config['sample_right']
        self.alpha = config['alpha']
        # self.pre_init_map_scheme = config['pre_init_map_scheme']
        # self.init_map_scheme = config['init_map_scheme']
        self.numIterMaxEmd = config['numIterMaxEmd']
        # self.n_hvgs = config['n_hvgs']
        # self.n_pcs = config['n_pcs']
        self.use_gpu = config['use_gpu']
        # self.sample_left_hvg_h5_save_path = config['sample_left_hvg_h5_save_path']
        # self.sample_right_hvg_h5_save_path = config['sample_right_hvg_h5_save_path']
        # self.data_folder_path = config['data_folder_path']
        # self.preprocessing = config['preprocessing']
        self.results_path = config['results_path']
        self.pi = config['pi']
        self.config_file_name = os.path.basename(self.config['config_path'])
        self.dissimilarity = config['dissimilarity']
        self.lambda_sinkhorn = config['lambda_sinkhorn']
        self.sinkhorn = config['sinkhorn']
        self.numInnerIterMax = config['numInnerIterMax']
        self.grid_search = config['grid_search']

        if config['adata_left_path'] != 'None':
            self.adata_left = sc.read(config['adata_left_path'])
            self.adata_right = sc.read(config['adata_right_path'])
        else:
            data_loader = DataLoader(config)

            dataset_map = data_loader.read_data(self.dataset)
            print('----------------------------------')
            print(self.sample_left)
            print(self.sample_right)
            print('----------------------------------')
            self.adata_left = dataset_map[self.sample_left]
            self.adata_right = dataset_map[self.sample_right]
        
        print("here:", self.adata_left.n_obs)

        # self.distribution = gamma
        # self.gamma_a = 0
        # self.gamma_loc = 0
        # self.gamma_scale = 0

        # if self.init_map_scheme == 'spatial':
        #     coor_1 = self.adata_left.obsm['spatial']

        #     with open(f'{self.results_path}/{self.dataset}/{self.config_file_name}/init_transform.json') as f:
        #         init_transform = json.load(f)
        #     if init_transform['mirror']:
        #         coor_1 = mirror(coor_1)
        #     coor_1 = rotate(coor_1, init_transform['init_rotation'])
        #     self.adata_left.obsm['spatial'] = coor_1
        
        config['cost_mat'] = np.load(config['cost_mat_path'])
        self.cost_mat = config['cost_mat']

        scale_coords(self.adata_left, key_name='spatial')
        scale_coords(self.adata_right, key_name='spatial')

        if config['QC']:
            QC(self.adata_left)
            QC(self.adata_right)

        self.fig_hist_rs, self.ax_hist_rs = plt.subplots()
        self.ax_hist_rs.set_xlabel('Remodeling score')
        self.ax_hist_rs.set_ylabel('Count')
        
    # def visualize_mapping(self):
    #     config_file_name = os.path.basename(self.config['config_path'])
    #     adata_left, adata_right, pi = self.adata_left, self.adata_right, self.pi
    #     os.makedirs(f'{self.results_path}/{self.dataset}/{config_file_name}/Mappings/', exist_ok=True)
    #     non_zero_count = np.count_nonzero(pi)
    #     plot_slice_pairwise_alignment_modified(adata_left, adata_right, pi, top=non_zero_count, save=True, save_dir=f'{self.results_path}/{self.dataset}/{config_file_name}/Mappings', name=f'{self.sample_left}_vs_{self.sample_right}_mapping.jpg')

    # def visualize_coeff_of_variation(self):
    #     f, ax = plt.subplots()
    #     f.set_size_inches(10, 10)
    #     points = ax.scatter(self.adata_right.obsm['spatial'][:, 0], self.adata_right.obsm['spatial'][:, 1], c=variation(self.pi, axis=0))
    #     f.colorbar(points)
    #     plt.close()

    def visualize_goodness_of_mapping(self, slice_pos='right', invert_x=False):
        if slice_pos == 'left': 
            adata = self.adata_left
            pi = self.pi
            cost_mat = self.cost_mat
            sample_name = self.sample_left
        else:
            adata = self.adata_right
            pi = self.pi.T
            cost_mat = self.cost_mat.T
            sample_name = self.sample_right

        score_mat = pi * cost_mat

        adata.obs['pathological_score'] = np.sum(score_mat, axis=1, dtype=np.float64) / (1 / adata.n_obs) * 100
        adata.obs['pathological_score'].to_csv(f'{self.results_path}/{self.dataset}/{self.config_file_name}/pathological_scores.csv')

        bins = 100
        plt.figure(figsize=(9, 9))
        plt.hist(adata.obs['pathological_score'].values, bins=bins)
        os.makedirs(f'{self.results_path}/{self.dataset}/{self.config_file_name}/Histograms/', exist_ok=True)
        plt.savefig(f'{self.results_path}/{self.dataset}/{self.config_file_name}/Histograms/{slice_pos}_pathological_score.jpg',format='jpg',dpi=350,bbox_inches='tight',pad_inches=0)
        plt.close()

        f, ax = plt.subplots()
        plt.figure(figsize=(9, 9))
        ax.axis('off')
        points = ax.scatter(-adata.obsm['spatial'][:, 0], -adata.obsm['spatial'][:, 1], s=10, c=-adata.obs['pathological_score'].values, cmap='plasma')
        if invert_x:
            f.gca().invert_xaxis()
        f.colorbar(points)
        config_file_name = os.path.basename(self.config['config_path'])
        os.makedirs(f'{self.results_path}/{self.dataset}/{config_file_name}/Goodness_score/', exist_ok=True)
        f.savefig(f'{self.results_path}/{self.dataset}/{config_file_name}/Goodness_score/{sample_name}_goodness_of_mapping.jpg',format='jpg',dpi=350,bbox_inches='tight',pad_inches=0)
        f.savefig(f'{self.results_path}/{self.dataset}/{config_file_name}/Goodness_score/{sample_name}_goodness_of_mapping.eps',format='eps',dpi=350,bbox_inches='tight',pad_inches=0)
        f.savefig(f'{self.results_path}/{self.dataset}/{config_file_name}/Goodness_score/{sample_name}_goodness_of_mapping.svg',format='svg',dpi=350,bbox_inches='tight',pad_inches=0)
        plt.close()

    # def get_p_value_wrt_gamma(self, val, a, loc, scale):
    #     cdf = gamma.cdf(val, a, loc=loc, scale=scale)
    #     p_value = 1 - cdf
    #     return p_value

    def divide_into_2_regions_wrt_goodness_score_and_find_DEG(self):
        # left_threshold = get_goodness_threshold_from_null_distribution(self.adata_left)
        # adata_to_be_synthesized = sc.read('../../../Data/King/Fixed_adatas/adata_Sham_1.h5ad')
        if self.config['adata_to_be_synthesized_path'] != 'None':
            adata_to_be_synthesized = sc.read(self.config['adata_to_be_synthesized_path'])
        else:
            adata_to_be_synthesized = self.adata_left.copy()

        adata_healthy_right = 'None'
        decompose = True
        if self.config['adata_healthy_right_path'] != 'None':
            decompose = False
            adata_healthy_right = sc.read(self.config['adata_healthy_right_path'])

        right_threshold = self.get_goodness_threshold_from_null_distribution(adata_to_be_synthesized, adata_healthy_right, decompose=decompose)

        # p_values = np.array(list(map(lambda x: self.get_p_value_wrt_gamma(x, self.gamma_a, self.gamma_loc, self.gamma_scale), self.adata_right.obs['pathological_score'].values)))
        # is_remodeled = np.array([0] * self.adata_right.n_obs)
        # is_remodeled[np.where(p_values < 0.05)[0]] = 1
        # np.save('../../../Workspace/SPaSE/local_data/is_remodeled.npy', is_remodeled)

        # plt.figure(figsize = (10, 10))
        # plt.axis('off')
        # plt.scatter(self.adata_right.obsm['spatial'][:, 0], self.adata_right.obsm['spatial'][:, 1], c=-is_remodeled, cmap='plasma')
        # plt.savefig(f'{self.results_path}/{self.dataset}/{self.config_file_name}/segmentation_based_of_fitted_gamma_distribution.jpg')

        print('Thresholds:', right_threshold)
        df_right_threshold = pd.DataFrame({'right_threshold': [right_threshold]})
        df_right_threshold.to_csv(f'{self.results_path}/{self.dataset}/{self.config_file_name}/thresholds.csv')

        self.adata_right.obs['region'] = 'bad'
        self.adata_right.obs.loc[self.adata_right.obs['pathological_score'] < right_threshold, 'region'] = 'good'
        self.adata_right.obs['region'] = self.adata_right.obs['region'].astype('category')

        plt.close('all')
        plt.figure(figsize = (10, 10))
        plt.axis('off')
        plt.scatter(self.adata_right.obsm['spatial'][:, 0], self.adata_right.obsm['spatial'][:, 1], c = list(map(lambda x: 1 if x=='good' else 0, pd.Categorical(self.adata_right.obs['region']))), cmap='plasma')
        plt.savefig(f'{self.results_path}/{self.dataset}/{self.config_file_name}/segmentation_based_on_discrete_distribution.jpg')
        plt.close()

        # os.makedirs(f'{self.results_path}/../local_data/{self.config_file_name}/Processed_adatas/', exist_ok=True)
        # self.adata_right.write(f'{self.results_path}/../local_data/{self.config_file_name}/Processed_adatas/adata_right_processed.h5ad')

        pi_right_to_left = self.pi.T
        region_col = self.adata_right.obs['region'].values
        idx_adata_right_bad = np.where(region_col == 'bad')[0]
        
        col_sum = pi_right_to_left[idx_adata_right_bad].sum(axis=0)
        mapped_bad_idx_left_int = np.where(col_sum != 0)[0]
        idx_barcodes = self.adata_left.obs.index[mapped_bad_idx_left_int]
        self.adata_left.obs['region_mapped'] = 'good'
        self.adata_left.obs.loc[idx_barcodes, 'region_mapped'] = 'bad'

        sns.histplot(self.adata_right.obs['pathological_score'].values, kde=True, color="blue", ax=self.ax_hist_rs, bins=100)
        self.ax_hist_rs.legend(['Left (H)', 'Right (D)'])

        self.fig_hist_rs.savefig(f'{self.results_path}/{self.dataset}/{self.config_file_name}/rs_distribution_both_both_samples.jpg')

        if self.grid_search:
            actual = self.adata_right.obs['is_remodeled_for_grid_search'].values
            predicted = np.array(list(map(lambda x: 1 if x == 'bad' else 0, self.adata_right.obs['region'].values)))
            # predicted_gamma = is_remodeled # prediction from the fitted gamma distribution

            F1_score = metrics.f1_score(actual, predicted)
            # F1_score_gamma = metrics.f1_score(actual, predicted_gamma)

            df_F1_score = pd.DataFrame({'F1_score': [F1_score]})
            df_F1_score.to_csv(f'{self.results_path}/{self.dataset}/{self.config_file_name}/F1_score.csv')

            # df_F1_score_gamma = pd.DataFrame({'F1_score': [F1_score_gamma]})
            # df_F1_score_gamma.to_csv(f'{self.results_path}/{self.dataset}/{self.config_file_name}/F1_score_gamma.csv')

        # os.makedirs(f'{self.results_path}/../local_data/{self.config_file_name}/Processed_adatas/', exist_ok=True)
        # self.adata_left.write(f'{self.results_path}/../local_data/{self.config_file_name}/Processed_adatas/adata_left_processed.h5ad')

    def get_goodness_threshold_from_null_distribution(self, adata, adata_2='None', decompose=True):
        print("\nSynthesizing the healthy sample\n")
        if decompose:
            adata_0, adata_1 = get_2hop_adatas(adata)
        else:
            adata_0 = adata
            adata_1 = adata_2

            print(adata_0.n_obs)
            print(adata_1.n_obs)

        backend = ot.backend.NumpyBackend()
        use_gpu = False
        if torch.cuda.is_available():
            backend = ot.backend.TorchBackend()
            use_gpu = True
        
        if decompose:
            cost_mat_path = f'{self.results_path}/../local_data/{self.dataset}/{self.sample_left}/cost_mat_{self.sample_left}_0_{self.sample_left}_1_{self.dissimilarity}.npy'
        else:
            cost_mat_path = f'{self.results_path}/../local_data/{self.dataset}/{self.sample_left}/cost_mat_Sham_1_Sham_2_{self.dissimilarity}.npy'
        os.makedirs(os.path.dirname(cost_mat_path), exist_ok=True)
        # pi_low_entropy_path = f'{self.results_path}/{self.dataset}/config_{self.dataset}_{self.sample_left}_vs_{self.sample_right}_js.json/Pis/{self.dataset}_synthetic_left_right_low_entropy.npy'
        plt.switch_backend('agg')
        if self.sinkhorn:
            pi_low_entropy = None
            pi, fgw_dist = paste_pairwise_align_modified(adata_0,
                                             adata_1,
                                             alpha=self.alpha,
                                             G_init=pi_low_entropy,
                                             numItermax=10000,
                                             dissimilarity=self.dissimilarity,
                                             sinkhorn=self.sinkhorn,
                                             lambda_sinkhorn=self.lambda_sinkhorn,
                                             cost_mat_path=cost_mat_path,
                                             return_obj=True,
                                             verbose=False,
                                             norm=True,
                                             backend=backend,
                                             use_gpu=use_gpu,
                                             numInnerItermax=self.numInnerIterMax)
            # np.save(f'{self.results_path}/../local_data/{self.dataset}/{self.sample_left}/pi_adata_0_vs_adata_1_alpha_{self.alpha}_lambda_{self.lambda_sinkhorn}.npy', pi)
            
        else:
            print('sinkhorn not used')
            pi = paste_pairwise_align_modified(adata_0, adata_1,alpha=self.alpha,G_init=None,numItermax=10000,dissimilarity=self.dissimilarity,sinkhorn=self.sinkhorn,cost_mat_path=cost_mat_path,verbose=False,backend=backend,use_gpu=use_gpu, norm=True, numItermaxEmd=self.numIterMaxEmd)

        cost_mat = np.load(cost_mat_path)


        distances_left, weights_left = compute_null_distribution(pi, cost_mat, 'left')
        # print('\n\ndistances_left', distances_left.min(), distances_left.max(), '\n\n')
        # a_left, loc_left, scale_left = self.distribution.fit(distances_left)
        # print(a_left, loc_left, scale_left)
        

        plt.figure(figsize=(9, 9))
        plt.tick_params(axis='x', labelsize=15)
        plt.tick_params(axis='y', labelsize=15)
        left_freqs = plt.hist(distances_left, bins=100)[0]
        
        os.makedirs(f'{self.results_path}/{self.dataset}/{self.config_file_name}/Histograms/', exist_ok=True)
        plt.savefig(f'{self.results_path}/{self.dataset}/{self.config_file_name}/Histograms/splitted_slice_left_pathological_score.jpg',format='jpg',dpi=350,bbox_inches='tight',pad_inches=0)

        distances_right, weights_right = compute_null_distribution(pi, cost_mat, 'right')
        print('\n\ndistances_right', distances_right.min(), distances_right.max(), '\n\n')

        # a_right, loc_right, scale_right = self.distribution.fit(distances_right)
        # print(a_right, loc_right, scale_right)

        plt.figure(figsize=(9, 9))
        plt.tick_params(axis='x', labelsize=15)
        plt.tick_params(axis='y', labelsize=15)
        right_freqs = plt.hist(distances_right, bins=100)[0]

        os.makedirs(f'{self.results_path}/{self.dataset}/{self.config_file_name}/Histograms/', exist_ok=True)
        plt.savefig(f'{self.results_path}/{self.dataset}/{self.config_file_name}/Histograms/splitted_slice_right_pathological_score.jpg',format='jpg',dpi=350,bbox_inches='tight',pad_inches=0)

        p_value = stats.kstest(left_freqs, right_freqs)[1]
        print('KS test pvalue:', p_value)

        distances_both = np.array(list(distances_left) + list(distances_right))
        weights_both = np.array(list(weights_left) + list(weights_right))
        
        # a_both, loc_both, scale_both = self.distribution.fit(distances_both)

        significance_threshold = 0.95

        # self.gamma_a = a_both
        # self.gamma_loc = loc_both
        # self.gamma_scale = scale_both
    
        bin_values_both = plt.hist(distances_both, weights=weights_both, bins = 100)
        pd.DataFrame({'Synthetic_spot_dist_both': distances_both}).to_csv(f'{self.results_path}/{self.dataset}/{self.config_file_name}/synthetic_spot_distances_both.csv')
        freqs = bin_values_both[0]
        print(max(bin_values_both[1]))
        sum_both = 0
        sum_tot = sum(freqs)
        for i in range(len(freqs)):
            sum_both += freqs[i]
            if sum_both/sum_tot > significance_threshold:
                both_threshold = bin_values_both[1][i]
                break
        sns.histplot(distances_both, kde=True, color="red", ax=self.ax_hist_rs, bins=100)
        
        return both_threshold