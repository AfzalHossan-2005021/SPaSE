import scanpy as sc
import seaborn as sns
import pandas as pd
import numpy as np
from anndata import AnnData
from scipy import sparse



class DataLoader():
    def __init__(self, config):
        self.data_folder_path = config['data_folder_path']

    def check_for_negative(self, adata):
        try:
            if adata.X.min() < 0:
                print('Found negative value is count matrix')
        except:
            if adata.X.toarray().min() < 0:
                print('Found negative value is count matrix')

    def read_data(self, dataset):
        if dataset == "King":
            adata_sham_1 = sc.read_visium(
                f'{self.data_folder_path}/King/V_sham_1_Np_2/', count_file='filtered_feature_bc_matrix.h5')
            adata_sham_1.var_names_make_unique()

            adata_1hr = sc.read_visium(
                f'{self.data_folder_path}/King/V_1hr/', count_file='filtered_feature_bc_matrix.h5')
            adata_1hr.var_names_make_unique()

            adata_4hr = sc.read_visium(
                f'{self.data_folder_path}/King/V_4hr/', count_file='filtered_feature_bc_matrix.h5')
            adata_4hr.var_names_make_unique()

            adata_d3_1 = sc.read_visium(
                f'{self.data_folder_path}/King/V_d3_1/', count_file='filtered_feature_bc_matrix.h5')
            adata_d3_1.var_names_make_unique()

            adata_d3_2 = sc.read_visium(
                f'{self.data_folder_path}/King/V_d3_2/', count_file='filtered_feature_bc_matrix.h5')
            adata_d3_2.var_names_make_unique()

            adata_d3_3 = sc.read_visium(
                f'{self.data_folder_path}/King/V_d3_3/', count_file='filtered_feature_bc_matrix.h5')
            adata_d3_3.var_names_make_unique()

            adata_d7_2 = sc.read_visium(
                f'{self.data_folder_path}/King/V_d7_2/', count_file='filtered_feature_bc_matrix.h5')
            adata_d7_2.var_names_make_unique()

            adata_d7_3 = sc.read_visium(
                f'{self.data_folder_path}/King/V_d7_3/', count_file='filtered_feature_bc_matrix.h5')
            adata_d7_3.var_names_make_unique()

            return {'Sham_1': adata_sham_1, '1hr': adata_1hr, '4hr': adata_4hr, 'D3_1': adata_d3_1, 'D3_3': adata_d3_3, 'D7_2': adata_d7_2, 'D7_3': adata_d7_3}

        elif dataset == "King_fixed_perturbed":
            # adata_sham_1 = sc.read_visium(
            #     f'{self.data_folder_path}/King/V_sham_1_Np_2/', count_file='filtered_feature_bc_matrix.h5')
            # adata_sham_1.var_names_make_unique()

            adata_sham_1 = sc.read(
                f'{self.data_folder_path}/King/Fixed_adatas/adata_Sham_1.h5ad')
            adata_sham_1.var_names_make_unique()
            # x = adata_sham_1.obsm['spatial'][:, 0]
            # y = adata_sham_1.obsm['spatial'][:, 1]
            # adata_sham_1.obsm['spatial'][:, 0] = x / x.max()
            # adata_sham_1.obsm['spatial'][:, 1] = y / y.max()

            # adata_perturbed_Sham_1_cov_scale_1_type_2 = sc.read(
            #     f'{self.data_folder_path}/King/Perturbed_adatas/adata_perturbed_Sham_1_cov_scale_1_type_2.h5ad')
            # adata_perturbed_Sham_1_cov_scale_1_type_2.var_names_make_unique()

            # adata_perturbed_Sham_1_cov_scale_10_type_2 = sc.read(
            #     f'{self.data_folder_path}/King/Perturbed_adatas/adata_perturbed_Sham_1_cov_scale_10_type_2.h5ad')
            # adata_perturbed_Sham_1_cov_scale_10_type_2.var_names_make_unique()

            # adata_perturbed_Sham_1_cov_scale_100_type_2 = sc.read(
            #     f'{self.data_folder_path}/King/Perturbed_adatas/adata_perturbed_Sham_1_cov_scale_100_type_2.h5ad')
            # adata_perturbed_Sham_1_cov_scale_100_type_2.var_names_make_unique()



            adata_perturbed_Sham_1_cov_scale_1_type_2_1 = sc.read(
                f'{self.data_folder_path}/King/Perturbed_adatas/adata_perturbed_Sham_1_cov_scale_1_type_2_1.h5ad')
            adata_perturbed_Sham_1_cov_scale_1_type_2_1.var_names_make_unique()

            adata_perturbed_Sham_1_cov_scale_1_type_2_2 = sc.read(
                f'{self.data_folder_path}/King/Perturbed_adatas/adata_perturbed_Sham_1_cov_scale_1_type_2_2.h5ad')
            adata_perturbed_Sham_1_cov_scale_1_type_2_2.var_names_make_unique()

            adata_perturbed_Sham_1_cov_scale_1_type_2_3 = sc.read(
                f'{self.data_folder_path}/King/Perturbed_adatas/adata_perturbed_Sham_1_cov_scale_1_type_2_3.h5ad')
            adata_perturbed_Sham_1_cov_scale_1_type_2_3.var_names_make_unique()

            adata_perturbed_Sham_1_cov_scale_1_type_2_4 = sc.read(
                f'{self.data_folder_path}/King/Perturbed_adatas/adata_perturbed_Sham_1_cov_scale_1_type_2_4.h5ad')
            adata_perturbed_Sham_1_cov_scale_1_type_2_4.var_names_make_unique()

            adata_perturbed_Sham_1_cov_scale_1_type_2_5 = sc.read(
                f'{self.data_folder_path}/King/Perturbed_adatas/adata_perturbed_Sham_1_cov_scale_1_type_2_5.h5ad')
            adata_perturbed_Sham_1_cov_scale_1_type_2_5.var_names_make_unique()

            # adata_perturbed_Sham_1_cov_scale_1_type_2_disk_test_1 = sc.read(
            #     f'{self.data_folder_path}/King/Perturbed_adatas/adata_perturbed_Sham_1_cov_scale_1_type_2_disk_test_1.h5ad')
            # adata_perturbed_Sham_1_cov_scale_1_type_2_disk_test_1.var_names_make_unique()




            adata_perturbed_Sham_1_cov_scale_1_type_3 = sc.read(
                f'{self.data_folder_path}/King/Perturbed_adatas/adata_perturbed_Sham_1_cov_scale_1_type_3.h5ad')
            adata_perturbed_Sham_1_cov_scale_1_type_3.var_names_make_unique()
            # x = adata_perturbed_Sham_1_cov_scale_1_type_3.obsm['spatial'][:, 0]
            # y = adata_perturbed_Sham_1_cov_scale_1_type_3.obsm['spatial'][:, 1]
            # adata_perturbed_Sham_1_cov_scale_1_type_3.obsm['spatial'][:, 0] = x / x.max()
            # adata_perturbed_Sham_1_cov_scale_1_type_3.obsm['spatial'][:, 1] = y / y.max()

            adata_perturbed_Sham_1_cov_scale_10_type_3 = sc.read(
                f'{self.data_folder_path}/King/Perturbed_adatas/adata_perturbed_Sham_1_cov_scale_10_type_3.h5ad')
            adata_perturbed_Sham_1_cov_scale_10_type_3.var_names_make_unique()

            adata_perturbed_Sham_1_cov_scale_100_type_3 = sc.read(
                f'{self.data_folder_path}/King/Perturbed_adatas/adata_perturbed_Sham_1_cov_scale_100_type_3.h5ad')
            adata_perturbed_Sham_1_cov_scale_100_type_3.var_names_make_unique()

            # adata_sham_1_perturbed = sc.read(f'{self.data_folder_path}/King/Perturbed_adatas/adata_perturbed_Sham_1_cov_scale_1_type_1.h5ad')
            # adata_sham_1_perturbed.var_names_make_unique()
            
            

            return {
                'Sham_1': adata_sham_1,
                'Sham_1_perturbed_cov_scale_1_type_2_1': adata_perturbed_Sham_1_cov_scale_1_type_2_1,
                'Sham_1_perturbed_cov_scale_1_type_2_2': adata_perturbed_Sham_1_cov_scale_1_type_2_2,
                'Sham_1_perturbed_cov_scale_1_type_2_3': adata_perturbed_Sham_1_cov_scale_1_type_2_3,
                'Sham_1_perturbed_cov_scale_1_type_2_4': adata_perturbed_Sham_1_cov_scale_1_type_2_4,
                'Sham_1_perturbed_cov_scale_1_type_2_5': adata_perturbed_Sham_1_cov_scale_1_type_2_5,
                # 'Sham_1_perturbed_cov_scale_1_type_2': adata_perturbed_Sham_1_cov_scale_1_type_2,
                # 'Sham_1_perturbed_cov_scale_10_type_2': adata_perturbed_Sham_1_cov_scale_10_type_2,
                # 'Sham_1_perturbed_cov_scale_100_type_2': adata_perturbed_Sham_1_cov_scale_100_type_2,
                'Sham_1_perturbed_cov_scale_1_type_3': adata_perturbed_Sham_1_cov_scale_1_type_3,
                'Sham_1_perturbed_cov_scale_10_type_3': adata_perturbed_Sham_1_cov_scale_10_type_3,
                'Sham_1_perturbed_cov_scale_100_type_3': adata_perturbed_Sham_1_cov_scale_100_type_3,
                # 'Sham_1_perturbed_1_cov_scale_1_type_2_disk_test_1': adata_perturbed_Sham_1_cov_scale_1_type_2_disk_test_1
                }
            # adata_d3_1 = sc.read(f'{self.data_folder_path}/King/Perturbed_adatas/adata_perturbed_D3_1.h5ad')
            # adata_d3_3 = sc.read(f'{self.data_folder_path}/King/Perturbed_adatas/adata_perturbed_D3_3.h5ad')
            # adata_d7_2 = sc.read(f'{self.data_folder_path}/King/Perturbed_adatas/adata_perturbed_D7_2.h5ad')
            # adata_d7_3 = sc.read(f'{self.data_folder_path}/King/Perturbed_adatas/adata_perturbed_D7_3.h5ad')

            # return {'Sham_1': adata_sham_1, 'Sham_1_perturbed': adata_sham_1_perturbed, 'D3_1': adata_d3_1, 'D3_3': adata_d3_3, 'D7_2': adata_d7_2, 'D7_3': adata_d7_3}
        elif dataset == "King_fixed_perturbed_initialized":
            adata_perturbed_Sham_1_cov_scale_1_type_2_disk_test_1 = sc.read(
                f'{self.data_folder_path}/King/Perturbed_adatas/adata_perturbed_Sham_1_cov_scale_1_type_2_disk_test_1.h5ad')
            adata_perturbed_Sham_1_cov_scale_1_type_2_disk_test_1.var_names_make_unique()

            adata_sham_1 = sc.read(
                f'{self.data_folder_path}/King/Fixed_adatas/adata_Sham_1.h5ad')
            adata_sham_1.var_names_make_unique()

            return {
                'Sham_1': adata_sham_1,
                'Sham_1_perturbed_1_cov_scale_1_type_2_disk_test_1': adata_perturbed_Sham_1_cov_scale_1_type_2_disk_test_1
            }

        elif dataset == 'King_perturbed_grid_search':

            adata_sham_1 = sc.read(
                f'{self.data_folder_path}/King/Fixed_adatas/adata_Sham_1.h5ad')
            adata_sham_1.var_names_make_unique()

            adata_sham_1_copy = sc.read(
                f'{self.data_folder_path}/King/Fixed_adatas/adata_Sham_1.h5ad')
            adata_sham_1_copy.var_names_make_unique()


            adata_perturbed_Sham_1_cov_scale_1_type_3_double_remodeling_1 = sc.read(
                f'{self.data_folder_path}/King/Perturbed_adatas/adata_perturbed_Sham_1_cov_scale_1_type_3_double_remodeling_1.h5ad')
            adata_perturbed_Sham_1_cov_scale_1_type_3_double_remodeling_1.var_names_make_unique()

            adata_perturbed_Sham_1_cov_scale_1_type_3_double_remodeling_2 = sc.read(
                f'{self.data_folder_path}/King/Perturbed_adatas/adata_perturbed_Sham_1_cov_scale_1_type_3_double_remodeling_2.h5ad')
            adata_perturbed_Sham_1_cov_scale_1_type_3_double_remodeling_2.var_names_make_unique()

            adata_perturbed_Sham_1_cov_scale_1_type_3_double_remodeling_3 = sc.read(
                f'{self.data_folder_path}/King/Perturbed_adatas/adata_perturbed_Sham_1_cov_scale_1_type_3_double_remodeling_3.h5ad')
            adata_perturbed_Sham_1_cov_scale_1_type_3_double_remodeling_3.var_names_make_unique()

            adata_perturbed_Sham_1_cov_scale_1_type_3_double_remodeling_4 = sc.read(
                f'{self.data_folder_path}/King/Perturbed_adatas/adata_perturbed_Sham_1_cov_scale_1_type_3_double_remodeling_4.h5ad')
            adata_perturbed_Sham_1_cov_scale_1_type_3_double_remodeling_4.var_names_make_unique()

            adata_perturbed_Sham_1_cov_scale_1_type_3_double_remodeling_5 = sc.read(
                f'{self.data_folder_path}/King/Perturbed_adatas/adata_perturbed_Sham_1_cov_scale_1_type_3_double_remodeling_5.h5ad')
            adata_perturbed_Sham_1_cov_scale_1_type_3_double_remodeling_5.var_names_make_unique()


            adata_perturbed_Sham_1_cov_scale_1_type_2_double_remodeling_all_non_zero_genes_1 = sc.read(
                f'{self.data_folder_path}/King/Perturbed_adatas/adata_perturbed_Sham_1_cov_scale_1_type_2_double_remodeling_all_non_zero_genes_1.h5ad')
            adata_perturbed_Sham_1_cov_scale_1_type_2_double_remodeling_all_non_zero_genes_1.var_names_make_unique()

            adata_perturbed_Sham_1_cov_scale_1_type_2_double_remodeling_all_non_zero_genes_2 = sc.read(
                f'{self.data_folder_path}/King/Perturbed_adatas/adata_perturbed_Sham_1_cov_scale_1_type_2_double_remodeling_all_non_zero_genes_2.h5ad')
            adata_perturbed_Sham_1_cov_scale_1_type_2_double_remodeling_all_non_zero_genes_2.var_names_make_unique()

            adata_perturbed_Sham_1_cov_scale_1_type_2_double_remodeling_all_non_zero_genes_3 = sc.read(
                f'{self.data_folder_path}/King/Perturbed_adatas/adata_perturbed_Sham_1_cov_scale_1_type_2_double_remodeling_all_non_zero_genes_3.h5ad')
            adata_perturbed_Sham_1_cov_scale_1_type_2_double_remodeling_all_non_zero_genes_3.var_names_make_unique()


            return {
                'Sham_1': adata_sham_1,
                'Sham_1_copy': adata_sham_1_copy,
                'Sham_1_perturbed_cov_scale_1_type_3_double_remodeling_1': adata_perturbed_Sham_1_cov_scale_1_type_3_double_remodeling_1,
                'Sham_1_perturbed_cov_scale_1_type_3_double_remodeling_2': adata_perturbed_Sham_1_cov_scale_1_type_3_double_remodeling_2,
                'Sham_1_perturbed_cov_scale_1_type_3_double_remodeling_3': adata_perturbed_Sham_1_cov_scale_1_type_3_double_remodeling_3,
                'Sham_1_perturbed_cov_scale_1_type_3_double_remodeling_4': adata_perturbed_Sham_1_cov_scale_1_type_3_double_remodeling_4,
                'Sham_1_perturbed_cov_scale_1_type_3_double_remodeling_5': adata_perturbed_Sham_1_cov_scale_1_type_3_double_remodeling_5,

                'Sham_1_perturbed_cov_scale_1_type_2_double_remodeling_all_non_zero_genes_1': adata_perturbed_Sham_1_cov_scale_1_type_2_double_remodeling_all_non_zero_genes_1,
                'Sham_1_perturbed_cov_scale_1_type_2_double_remodeling_all_non_zero_genes_2': adata_perturbed_Sham_1_cov_scale_1_type_2_double_remodeling_all_non_zero_genes_2,
                'Sham_1_perturbed_cov_scale_1_type_2_double_remodeling_all_non_zero_genes_3': adata_perturbed_Sham_1_cov_scale_1_type_2_double_remodeling_all_non_zero_genes_3,
                }

        elif dataset == 'King_fixed':
            adata_1hr = sc.read(
                f'{self.data_folder_path}/King/Fixed_adatas/adata_1hr.h5ad')
            adata_1hr.var_names_make_unique()

            adata_4hr = sc.read(
                f'{self.data_folder_path}/King/Fixed_adatas/adata_4hr.h5ad')
            adata_4hr.var_names_make_unique()

            adata_D3_1 = sc.read(
                f'{self.data_folder_path}/King/Fixed_adatas/adata_D3_1.h5ad')
            adata_D3_1.var_names_make_unique()

            adata_D3_3 = sc.read(
                f'{self.data_folder_path}/King/Fixed_adatas/adata_D3_3.h5ad')
            adata_D3_3.var_names_make_unique()

            adata_D7_2 = sc.read(
                f'{self.data_folder_path}/King/Fixed_adatas/adata_D7_2.h5ad')
            adata_D7_2.var_names_make_unique()

            adata_D7_3 = sc.read(
                f'{self.data_folder_path}/King/Fixed_adatas/adata_D7_3.h5ad')
            adata_D7_3.var_names_make_unique()

            adata_Sham_1 = sc.read(
                f'{self.data_folder_path}/King/Fixed_adatas/adata_Sham_1.h5ad')
            adata_Sham_1.var_names_make_unique()

            return {
                'Sham_1': adata_Sham_1, '1hr': adata_1hr, '4hr': adata_4hr, 'D3_1': adata_D3_1, 'D3_3': adata_D3_3, 'D7_2': adata_D7_2, 'D7_3': adata_D7_3
            }
        elif dataset == 'Human_heart':
            adata_control_P1 = sc.read(
                f'{self.data_folder_path}/Human_heart/control_P1.h5ad')
            adata_control_P1.var_names_make_unique()
            adata_control_P1.obsm['spatial'] = adata_control_P1.obsm['X_spatial']

            adata_control_P7 = sc.read(
                f'{self.data_folder_path}/Human_heart/control_P7.h5ad')
            adata_control_P7.var_names_make_unique()
            adata_control_P7.obsm['spatial'] = adata_control_P7.obsm['X_spatial']

            adata_control_P8 = sc.read(
                f'{self.data_folder_path}/Human_heart/control_P8.h5ad')
            adata_control_P8.var_names_make_unique()
            adata_control_P8.obsm['spatial'] = adata_control_P8.obsm['X_spatial']

            adata_control_P17 = sc.read(
                f'{self.data_folder_path}/Human_heart/control_P17.h5ad')
            adata_control_P17.var_names_make_unique()
            adata_control_P17.obsm['spatial'] = adata_control_P17.obsm['X_spatial']

            adata_IZ_P3 = sc.read(
                f'{self.data_folder_path}/Human_heart/IZ_P3.h5ad')
            adata_IZ_P3.var_names_make_unique()
            adata_IZ_P3.obsm['spatial'] = adata_IZ_P3.obsm['X_spatial']
            
            adata_RZ_BZ_P3 = sc.read(
                f'{self.data_folder_path}/Human_heart/RZ_BZ_P3.h5ad')
            adata_RZ_BZ_P3.var_names_make_unique()
            adata_RZ_BZ_P3.obsm['spatial'] = adata_RZ_BZ_P3.obsm['X_spatial']
            
            adata_RZ_P9 = sc.read(
                f'{self.data_folder_path}/Human_heart/RZ_P9.h5ad')
            adata_RZ_P9.var_names_make_unique()
            adata_RZ_P9.obsm['spatial'] = adata_RZ_P9.obsm['X_spatial']

            

            return {
                'adata_control_P1': adata_control_P1,
                'adata_control_P7': adata_control_P7,
                'adata_control_P8': adata_control_P8,
                'adata_control_P17': adata_control_P17,
                'adata_IZ_P3': adata_IZ_P3,
                'adata_RZ_BZ_P3': adata_RZ_BZ_P3,
                'adata_RZ_P9': adata_RZ_P9,
            }
        elif dataset == 'King_splitted':
            adata_0 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Workspace/PASTE_modified/local_data/King_pertubed_grid_search/Splitted_adatas/adata_0.h5ad')
            adata_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Workspace/PASTE_modified/local_data/King_pertubed_grid_search/Splitted_adatas/adata_1.h5ad')
            return {
                'adata_0': adata_0,
                'adata_1': adata_1
            }
        elif dataset == 'King_perturbed_mvn_with_disk':
            adata_Sham_1 = sc.read(
                f'{self.data_folder_path}/King/Fixed_adatas/adata_Sham_1.h5ad')
            adata_Sham_1.var_names_make_unique()
            data_path = '/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas'

            adata_Sham_1_perturbed_mvn_with_disk_1 = sc.read(f'{data_path}/adata_perturbed_Sham_1_perturbed_mvn_cov_scale_1_type_2_double_remodeling_all_non_zero_genes_1.h5ad')
            adata_Sham_1_perturbed_mvn_with_disk_1.var_names_make_unique()
            adata_Sham_1_perturbed_mvn_with_disk_2 = sc.read(f'{data_path}/adata_perturbed_Sham_1_perturbed_mvn_cov_scale_1_type_2_double_remodeling_all_non_zero_genes_2.h5ad')
            adata_Sham_1_perturbed_mvn_with_disk_2.var_names_make_unique()
            adata_Sham_1_perturbed_mvn_with_disk_3 = sc.read(f'{data_path}/adata_perturbed_Sham_1_perturbed_mvn_cov_scale_1_type_2_double_remodeling_all_non_zero_genes_3.h5ad')
            adata_Sham_1_perturbed_mvn_with_disk_3.var_names_make_unique()

            return {
                'Sham_1': adata_Sham_1,
                "Sham_1_mvn_1": adata_Sham_1_perturbed_mvn_with_disk_1,
                "Sham_1_mvn_2": adata_Sham_1_perturbed_mvn_with_disk_2,
                "Sham_1_mvn_3": adata_Sham_1_perturbed_mvn_with_disk_3,
            }
        elif dataset == 'King_perturbed_mvn_all_gene_with_disk':
            adata_Sham_1 = sc.read(
                f'{self.data_folder_path}/King/Fixed_adatas/adata_Sham_1.h5ad')
            adata_Sham_1.var_names_make_unique()
            data_path = '/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas'

            adata_Sham_1_perturbed_mvn_all_gene_with_disk_1 = sc.read(f'{data_path}/adata_perturbed_Sham_1_perturbed_mvn_all_gene_cov_scale_1_type_2_double_remodeling_all_non_zero_genes_1.h5ad')
            adata_Sham_1_perturbed_mvn_all_gene_with_disk_1.var_names_make_unique()
            adata_Sham_1_perturbed_mvn_all_gene_with_disk_2 = sc.read(f'{data_path}/adata_perturbed_Sham_1_perturbed_mvn_all_gene_cov_scale_1_type_2_double_remodeling_all_non_zero_genes_2.h5ad')
            adata_Sham_1_perturbed_mvn_all_gene_with_disk_2.var_names_make_unique()
            adata_Sham_1_perturbed_mvn_all_gene_with_disk_3 = sc.read(f'{data_path}/adata_perturbed_Sham_1_perturbed_mvn_all_gene_cov_scale_1_type_2_double_remodeling_all_non_zero_genes_3.h5ad')
            adata_Sham_1_perturbed_mvn_all_gene_with_disk_3.var_names_make_unique()

            # adata_Sham_1_perturbed_mvn = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/All_spot_multivariate/adata_Sham_1.h5ad')
            # adata_Sham_1_perturbed_mvn.var_names_make_unique()

            return {
                'Sham_1': adata_Sham_1,
                "Sham_1_mvn_all_gene_1": adata_Sham_1_perturbed_mvn_all_gene_with_disk_1,
                "Sham_1_mvn_all_gene_2": adata_Sham_1_perturbed_mvn_all_gene_with_disk_2,
                "Sham_1_mvn_all_gene_3": adata_Sham_1_perturbed_mvn_all_gene_with_disk_3
                # 'Sham_1_perturbed_mvn': adata_Sham_1_perturbed_mvn,
            }
        elif dataset == 'King_perturbed_mvn_scale_10_all_gene_with_disk':
            adata_Sham_1 = sc.read(
                f'{self.data_folder_path}/King/Fixed_adatas/adata_Sham_1.h5ad')
            adata_Sham_1.var_names_make_unique()
            data_path = '/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas'

            adata_Sham_1_perturbed_mvn_all_gene_with_disk_1 = sc.read(f'{data_path}/adata_perturbed_Sham_1_perturbed_mvn_all_gene_cov_scale_10_type_2_double_remodeling_all_non_zero_genes_1.h5ad')
            adata_Sham_1_perturbed_mvn_all_gene_with_disk_1.var_names_make_unique()
            adata_Sham_1_perturbed_mvn_all_gene_with_disk_2 = sc.read(f'{data_path}/adata_perturbed_Sham_1_perturbed_mvn_all_gene_cov_scale_10_type_2_double_remodeling_all_non_zero_genes_2.h5ad')
            adata_Sham_1_perturbed_mvn_all_gene_with_disk_2.var_names_make_unique()
            adata_Sham_1_perturbed_mvn_all_gene_with_disk_3 = sc.read(f'{data_path}/adata_perturbed_Sham_1_perturbed_mvn_all_gene_cov_scale_10_type_2_double_remodeling_all_non_zero_genes_3.h5ad')
            adata_Sham_1_perturbed_mvn_all_gene_with_disk_3.var_names_make_unique()

            # adata_Sham_1_perturbed_mvn = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/All_spot_multivariate/adata_Sham_1.h5ad')
            # adata_Sham_1_perturbed_mvn.var_names_make_unique()

            return {
                'Sham_1': adata_Sham_1,
                "Sham_1_mvn_all_gene_1": adata_Sham_1_perturbed_mvn_all_gene_with_disk_1,
                "Sham_1_mvn_all_gene_2": adata_Sham_1_perturbed_mvn_all_gene_with_disk_2,
                "Sham_1_mvn_all_gene_3": adata_Sham_1_perturbed_mvn_all_gene_with_disk_3
                # 'Sham_1_perturbed_mvn': adata_Sham_1_perturbed_mvn,
            }
        elif dataset == 'King_perturbed_scale_10_non_zero_gene_with_disk':
            adata_Sham_1 = sc.read(
                f'{self.data_folder_path}/King/Fixed_adatas/adata_Sham_1.h5ad')
            adata_Sham_1.var_names_make_unique()
            data_path = '/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas'

            adata_1 = sc.read(f'{data_path}/adata_perturbed_Sham_1_cov_scale_10_type_2_double_remodeling_all_non_zero_genes_1.h5ad')
            adata_1.var_names_make_unique()

            adata_2 = sc.read(f'{data_path}/adata_perturbed_Sham_1_cov_scale_10_type_2_double_remodeling_all_non_zero_genes_2.h5ad')
            adata_2.var_names_make_unique()

            adata_3 = sc.read(f'{data_path}/adata_perturbed_Sham_1_cov_scale_10_type_2_double_remodeling_all_non_zero_genes_3.h5ad')
            adata_3.var_names_make_unique()

            # adata_Sham_1_perturbed_mvn = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/All_spot_multivariate/adata_Sham_1.h5ad')
            # adata_Sham_1_perturbed_mvn.var_names_make_unique()

            return {
                'Sham_1': adata_Sham_1,
                "Sham_1_scale_10_non_zero_genes_1": adata_1,
                "Sham_1_scale_10_non_zero_genes_2": adata_2,
                "Sham_1_scale_10_non_zero_genes_3": adata_3,
                # 'Sham_1_perturbed_mvn': adata_Sham_1_perturbed_mvn,
            }
        elif dataset == 'King_mvn_vs_mvn_with_disk':
            adata_mvn = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_Sham_1_perturbed_mvn_2.h5ad')
            self.check_for_negative(adata_mvn)

            adata_mvn_with_disk_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_Sham_1_perturbed_mvn_3_cov_scale_1_type_2_double_remodeling_hvg_1.h5ad')
            self.check_for_negative(adata_mvn_with_disk_1)
            adata_mvn_with_disk_2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_Sham_1_perturbed_mvn_3_cov_scale_1_type_2_double_remodeling_hvg_2.h5ad')
            self.check_for_negative(adata_mvn_with_disk_2)
            adata_mvn_with_disk_3 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_Sham_1_perturbed_mvn_3_cov_scale_1_type_2_double_remodeling_hvg_3.h5ad')
            self.check_for_negative(adata_mvn_with_disk_3)

            adata_mvn_with_disk_1_scale_10 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_Sham_1_perturbed_mvn_3_cov_scale_10_type_2_double_remodeling_hvg_1.h5ad')
            self.check_for_negative(adata_mvn_with_disk_1_scale_10)
            adata_mvn_with_disk_2_scale_10 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_Sham_1_perturbed_mvn_3_cov_scale_10_type_2_double_remodeling_hvg_2.h5ad')
            self.check_for_negative(adata_mvn_with_disk_2_scale_10)
            adata_mvn_with_disk_3_scale_10 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_Sham_1_perturbed_mvn_3_cov_scale_10_type_2_double_remodeling_hvg_3.h5ad')
            self.check_for_negative(adata_mvn_with_disk_3_scale_10)

            return {
                'Sham_1_mvn': adata_mvn,
                'Sham_1_mvn_disk_perturbed_1': adata_mvn_with_disk_1,
                'Sham_1_mvn_disk_perturbed_2': adata_mvn_with_disk_2,
                'Sham_1_mvn_disk_perturbed_3': adata_mvn_with_disk_3,
                'Sham_1_mvn_disk_perturbed_1_scale_10': adata_mvn_with_disk_1_scale_10,
                'Sham_1_mvn_disk_perturbed_2_scale_10': adata_mvn_with_disk_2_scale_10,
                'Sham_1_mvn_disk_perturbed_3_scale_10': adata_mvn_with_disk_3_scale_10,
            }
        elif dataset == 'King_mvn_vs_mvn_with_disk_scale_10':
            adata_mvn = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_Sham_1_perturbed_mvn_2.h5ad')
            self.check_for_negative(adata_mvn)

            # adata_mvn_with_disk_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_Sham_1_perturbed_mvn_3_cov_scale_1_type_2_double_remodeling_hvg_1.h5ad')
            # self.check_for_negative(adata_mvn_with_disk_1)
            # adata_mvn_with_disk_2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_Sham_1_perturbed_mvn_3_cov_scale_1_type_2_double_remodeling_hvg_2.h5ad')
            # self.check_for_negative(adata_mvn_with_disk_2)
            # adata_mvn_with_disk_3 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_Sham_1_perturbed_mvn_3_cov_scale_1_type_2_double_remodeling_hvg_3.h5ad')
            # self.check_for_negative(adata_mvn_with_disk_3)

            adata_mvn_with_disk_1_scale_10 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_Sham_1_perturbed_mvn_3_cov_scale_10_type_2_double_remodeling_hvg_1.h5ad')
            self.check_for_negative(adata_mvn_with_disk_1_scale_10)
            adata_mvn_with_disk_2_scale_10 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_Sham_1_perturbed_mvn_3_cov_scale_10_type_2_double_remodeling_hvg_2.h5ad')
            self.check_for_negative(adata_mvn_with_disk_2_scale_10)
            adata_mvn_with_disk_3_scale_10 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_Sham_1_perturbed_mvn_3_cov_scale_10_type_2_double_remodeling_hvg_3.h5ad')
            self.check_for_negative(adata_mvn_with_disk_3_scale_10)

            return {
                'Sham_1_mvn': adata_mvn,
                # 'Sham_1_mvn_disk_perturbed_1': adata_mvn_with_disk_1,
                # 'Sham_1_mvn_disk_perturbed_2': adata_mvn_with_disk_2,
                # 'Sham_1_mvn_disk_perturbed_3': adata_mvn_with_disk_3,
                'Sham_1_mvn_disk_perturbed_1_scale_10': adata_mvn_with_disk_1_scale_10,
                'Sham_1_mvn_disk_perturbed_2_scale_10': adata_mvn_with_disk_2_scale_10,
                'Sham_1_mvn_disk_perturbed_3_scale_10': adata_mvn_with_disk_3_scale_10,
            }
        elif dataset == 'King_synthetic_vs_synthetic_scale_1':
            adata_synthetic_0 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_synthetic_0.h5ad')

            adata_synthetic_1_perturbed_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_1_type_2_double_remodeling_synthetic_1_1.h5ad')
            adata_synthetic_1_perturbed_2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_1_type_2_double_remodeling_synthetic_1_2.h5ad')
            adata_synthetic_1_perturbed_3 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_1_type_2_double_remodeling_synthetic_1_3.h5ad')

            return {
                'Synthetic_0': adata_synthetic_0,
                'Synthetic_perturbed_1': adata_synthetic_1_perturbed_1,
                'Synthetic_perturbed_2': adata_synthetic_1_perturbed_2,
                'Synthetic_perturbed_3': adata_synthetic_1_perturbed_3,
            }
        elif dataset == 'King_synthetic_vs_synthetic_scale_2':
            adata_synthetic_0 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_synthetic_0.h5ad')

            adata_synthetic_1_perturbed_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_2_type_2_double_remodeling_synthetic_1_1.h5ad')
            adata_synthetic_1_perturbed_2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_2_type_2_double_remodeling_synthetic_1_2.h5ad')
            adata_synthetic_1_perturbed_3 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_2_type_2_double_remodeling_synthetic_1_3.h5ad')

            return {
                'Synthetic_0': adata_synthetic_0,
                'Synthetic_perturbed_1': adata_synthetic_1_perturbed_1,
                'Synthetic_perturbed_2': adata_synthetic_1_perturbed_2,
                'Synthetic_perturbed_3': adata_synthetic_1_perturbed_3,
            }
        elif dataset == 'King_synthetic_vs_synthetic_scale_3':
            adata_synthetic_0 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_synthetic_0.h5ad')

            adata_synthetic_1_perturbed_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_3_type_2_double_remodeling_synthetic_1_1.h5ad')
            adata_synthetic_1_perturbed_2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_3_type_2_double_remodeling_synthetic_1_2.h5ad')
            adata_synthetic_1_perturbed_3 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_3_type_2_double_remodeling_synthetic_1_3.h5ad')

            return {
                'Synthetic_0': adata_synthetic_0,
                'Synthetic_perturbed_1': adata_synthetic_1_perturbed_1,
                'Synthetic_perturbed_2': adata_synthetic_1_perturbed_2,
                'Synthetic_perturbed_3': adata_synthetic_1_perturbed_3,
            }
        elif dataset == 'King_synthetic_vs_synthetic_scale_4':
            adata_synthetic_0 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_synthetic_0.h5ad')

            adata_synthetic_1_perturbed_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_4_type_2_double_remodeling_synthetic_1_1.h5ad')
            adata_synthetic_1_perturbed_2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_4_type_2_double_remodeling_synthetic_1_2.h5ad')
            adata_synthetic_1_perturbed_3 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_4_type_2_double_remodeling_synthetic_1_3.h5ad')

            return {
                'Synthetic_0': adata_synthetic_0,
                'Synthetic_perturbed_1': adata_synthetic_1_perturbed_1,
                'Synthetic_perturbed_2': adata_synthetic_1_perturbed_2,
                'Synthetic_perturbed_3': adata_synthetic_1_perturbed_3,
            }
        elif dataset == 'King_synthetic_vs_synthetic_scale_5':
            adata_synthetic_0 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_synthetic_0.h5ad')

            adata_synthetic_1_perturbed_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_5_type_2_double_remodeling_synthetic_1_1.h5ad')
            adata_synthetic_1_perturbed_2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_5_type_2_double_remodeling_synthetic_1_2.h5ad')
            adata_synthetic_1_perturbed_3 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_5_type_2_double_remodeling_synthetic_1_3.h5ad')

            return {
                'Synthetic_0': adata_synthetic_0,
                'Synthetic_perturbed_1': adata_synthetic_1_perturbed_1,
                'Synthetic_perturbed_2': adata_synthetic_1_perturbed_2,
                'Synthetic_perturbed_3': adata_synthetic_1_perturbed_3,
            }
        elif dataset == 'King_synthetic_vs_synthetic_scale_6':
            adata_synthetic_0 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_synthetic_0.h5ad')

            adata_synthetic_1_perturbed_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_6_type_2_double_remodeling_synthetic_1_1.h5ad')
            adata_synthetic_1_perturbed_2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_6_type_2_double_remodeling_synthetic_1_2.h5ad')
            adata_synthetic_1_perturbed_3 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_6_type_2_double_remodeling_synthetic_1_3.h5ad')

            return {
                'Synthetic_0': adata_synthetic_0,
                'Synthetic_perturbed_1': adata_synthetic_1_perturbed_1,
                'Synthetic_perturbed_2': adata_synthetic_1_perturbed_2,
                'Synthetic_perturbed_3': adata_synthetic_1_perturbed_3,
            }
        elif dataset == 'King_synthetic_vs_synthetic_scale_7':
            adata_synthetic_0 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_synthetic_0.h5ad')

            adata_synthetic_1_perturbed_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_7_type_2_double_remodeling_synthetic_1_1.h5ad')
            adata_synthetic_1_perturbed_2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_7_type_2_double_remodeling_synthetic_1_2.h5ad')
            adata_synthetic_1_perturbed_3 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_7_type_2_double_remodeling_synthetic_1_3.h5ad')

            return {
                'Synthetic_0': adata_synthetic_0,
                'Synthetic_perturbed_1': adata_synthetic_1_perturbed_1,
                'Synthetic_perturbed_2': adata_synthetic_1_perturbed_2,
                'Synthetic_perturbed_3': adata_synthetic_1_perturbed_3,
            }
        elif dataset == 'King_synthetic_vs_synthetic_scale_8':
            adata_synthetic_0 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_synthetic_0.h5ad')

            adata_synthetic_1_perturbed_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_8_type_2_double_remodeling_synthetic_1_1.h5ad')
            adata_synthetic_1_perturbed_2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_8_type_2_double_remodeling_synthetic_1_2.h5ad')
            adata_synthetic_1_perturbed_3 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_8_type_2_double_remodeling_synthetic_1_3.h5ad')

            return {
                'Synthetic_0': adata_synthetic_0,
                'Synthetic_perturbed_1': adata_synthetic_1_perturbed_1,
                'Synthetic_perturbed_2': adata_synthetic_1_perturbed_2,
                'Synthetic_perturbed_3': adata_synthetic_1_perturbed_3,
            }
        elif dataset == 'King_synthetic_vs_synthetic_scale_9':
            adata_synthetic_0 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_synthetic_0.h5ad')

            adata_synthetic_1_perturbed_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_9_type_2_double_remodeling_synthetic_1_1.h5ad')
            adata_synthetic_1_perturbed_2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_9_type_2_double_remodeling_synthetic_1_2.h5ad')
            adata_synthetic_1_perturbed_3 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_9_type_2_double_remodeling_synthetic_1_3.h5ad')

            return {
                'Synthetic_0': adata_synthetic_0,
                'Synthetic_perturbed_1': adata_synthetic_1_perturbed_1,
                'Synthetic_perturbed_2': adata_synthetic_1_perturbed_2,
                'Synthetic_perturbed_3': adata_synthetic_1_perturbed_3,
            }
        elif dataset == 'King_synthetic_vs_synthetic_scale_10':
            adata_synthetic_0 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_synthetic_0.h5ad')

            adata_synthetic_1_perturbed_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_10_type_2_double_remodeling_synthetic_1_1.h5ad')
            adata_synthetic_1_perturbed_2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_10_type_2_double_remodeling_synthetic_1_2.h5ad')
            adata_synthetic_1_perturbed_3 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_cov_scale_10_type_2_double_remodeling_synthetic_1_3.h5ad')

            return {
                'Synthetic_0': adata_synthetic_0,
                'Synthetic_perturbed_1': adata_synthetic_1_perturbed_1,
                'Synthetic_perturbed_2': adata_synthetic_1_perturbed_2,
                'Synthetic_perturbed_3': adata_synthetic_1_perturbed_3,
            }
        elif dataset == 'King_synthetic_vs_synthetic_type_4':
            adata_synthetic_0 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_synthetic_0.h5ad')
            adata_synthetic_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_synthetic_1.h5ad')
            adata_synthetic_1_perturbed_at_827 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_synthetic_1_type_4_orig_perturb_centering_827.h5ad')

            adata_synthetic_1_perturbed_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_type_4_orig_perturb_centering_827_type_4_remodeling_std_mul_2_sample_num_1.h5ad')
            adata_synthetic_1_perturbed_2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_type_4_orig_perturb_centering_827_type_4_remodeling_std_mul_2_sample_num_2.h5ad')
            adata_synthetic_1_perturbed_3 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_type_4_orig_perturb_centering_827_type_4_remodeling_std_mul_2_sample_num_3.h5ad')

            adata_synthetic_1_scale_1_perturbed_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_type_4_orig_perturb_centering_827_type_4_remodeling_std_mul_1_sample_num_1.h5ad')
            adata_synthetic_1_scale_1_perturbed_2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_type_4_orig_perturb_centering_827_type_4_remodeling_std_mul_1_sample_num_2.h5ad')
            adata_synthetic_1_scale_1_perturbed_3 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_type_4_orig_perturb_centering_827_type_4_remodeling_std_mul_1_sample_num_3.h5ad')

            adata_synthetic_1_scale_0_1_perturbed_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_type_4_orig_perturb_centering_827_type_4_remodeling_std_mul_1_sample_num_1.h5ad')
            adata_synthetic_1_scale_0_1_perturbed_2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_type_4_orig_perturb_centering_827_type_4_remodeling_std_mul_1_sample_num_2.h5ad')
            adata_synthetic_1_scale_0_1_perturbed_3 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_synthetic_1_type_4_orig_perturb_centering_827_type_4_remodeling_std_mul_1_sample_num_3.h5ad')

            return {
                'Synthetic_0': adata_synthetic_0,
                'Synthetic_1': adata_synthetic_1,
                'Synthetic_1_orig_perturbation': adata_synthetic_1_perturbed_at_827,
                'Synthetic_perturbed_1': adata_synthetic_1_perturbed_1,
                'Synthetic_perturbed_2': adata_synthetic_1_perturbed_2,
                'Synthetic_perturbed_3': adata_synthetic_1_perturbed_3,
                'Synthetic_perturbed_scale_1_1': adata_synthetic_1_scale_1_perturbed_1,
                'Synthetic_perturbed_scale_1_2': adata_synthetic_1_scale_1_perturbed_2,
                'Synthetic_perturbed_scale_1_3': adata_synthetic_1_scale_1_perturbed_3,
                'Synthetic_perturbed_scale_0.1_1': adata_synthetic_1_scale_0_1_perturbed_1,
                'Synthetic_perturbed_scale_0.1_2': adata_synthetic_1_scale_0_1_perturbed_2,
                'Synthetic_perturbed_scale_0.1_3': adata_synthetic_1_scale_0_1_perturbed_3,
            }
        elif dataset == 'King_diseased_data_perturbed_type_4':
            adata_Sham_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_Sham_1.h5ad')

            adata_D3_1_rem = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D3_1_type_4_remodeling.h5ad')
            adata_D3_1_rem.var_names_make_unique()

            adata_D3_3_rem = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D3_3_type_4_remodeling.h5ad')
            adata_D3_3_rem.var_names_make_unique()

            adata_D7_2_rem = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D7_2_type_4_remodeling.h5ad')
            adata_D7_2_rem.var_names_make_unique()

            adata_D7_3_rem = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D7_3_type_4_remodeling.h5ad')
            adata_D7_3_rem.var_names_make_unique()

            return {
                'Sham_1': adata_Sham_1,
                'D3_1_rem': adata_D3_1_rem,
                'D3_3_rem': adata_D3_3_rem,
                'D7_2_rem': adata_D7_2_rem,
                'D7_3_rem': adata_D7_3_rem,
            }
        elif dataset == 'King_diseased_data_perturbed_std_mul_200_type_4':
            adata_Sham_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_Sham_1.h5ad')

            adata_D3_1_rem_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D3_1_type_4_remodeling_std_mul_200_sample_num_1.h5ad')
            adata_D3_1_rem_1.var_names_make_unique()
            adata_D3_1_rem_2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D3_1_type_4_remodeling_std_mul_200_sample_num_2.h5ad')
            adata_D3_1_rem_2.var_names_make_unique()
            adata_D3_1_rem_3 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D3_1_type_4_remodeling_std_mul_200_sample_num_3.h5ad')
            adata_D3_1_rem_3.var_names_make_unique()

            adata_D3_3_rem_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D3_3_type_4_remodeling_std_mul_200_sample_num_1.h5ad')
            adata_D3_3_rem_1.var_names_make_unique()
            adata_D3_3_rem_2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D3_3_type_4_remodeling_std_mul_200_sample_num_2.h5ad')
            adata_D3_3_rem_2.var_names_make_unique()
            adata_D3_3_rem_3 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D3_3_type_4_remodeling_std_mul_200_sample_num_3.h5ad')
            adata_D3_3_rem_3.var_names_make_unique()

            adata_D7_2_rem_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D7_2_type_4_remodeling_std_mul_200_sample_num_1.h5ad')
            adata_D7_2_rem_1.var_names_make_unique()
            adata_D7_2_rem_2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D7_2_type_4_remodeling_std_mul_200_sample_num_2.h5ad')
            adata_D7_2_rem_2.var_names_make_unique()
            adata_D7_2_rem_3 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D7_2_type_4_remodeling_std_mul_200_sample_num_3.h5ad')
            adata_D7_2_rem_3.var_names_make_unique()

            adata_D7_3_rem_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D7_3_type_4_remodeling_std_mul_200_sample_num_1.h5ad')
            adata_D7_3_rem_1.var_names_make_unique()
            adata_D7_3_rem_2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D7_3_type_4_remodeling_std_mul_200_sample_num_2.h5ad')
            adata_D7_3_rem_2.var_names_make_unique()
            adata_D7_3_rem_3 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D7_3_type_4_remodeling_std_mul_200_sample_num_3.h5ad')
            adata_D7_3_rem_3.var_names_make_unique()

            return {
                'Sham_1': adata_Sham_1,
                'D3_1_rem_1': adata_D3_1_rem_1,
                'D3_1_rem_2': adata_D3_1_rem_2,
                'D3_1_rem_3': adata_D3_1_rem_3,

                'D3_3_rem_1': adata_D3_3_rem_1,
                'D3_3_rem_2': adata_D3_3_rem_2,
                'D3_3_rem_3': adata_D3_3_rem_3,

                'D7_2_rem_1': adata_D7_2_rem_1,
                'D7_2_rem_2': adata_D7_2_rem_2,
                'D7_2_rem_3': adata_D7_2_rem_3,

                'D7_3_rem_1': adata_D7_3_rem_1,
                'D7_3_rem_2': adata_D7_3_rem_2,
                'D7_3_rem_3': adata_D7_3_rem_3,
            }
        elif dataset == 'King_diseased_data_perturbed_std_mul_2_type_4':
            adata_Sham_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_Sham_1.h5ad')

            adata_1hr_rem_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_1hr_type_4_remodeling_std_mul_2_sample_num_1.h5ad')
            adata_1hr_rem_1.var_names_make_unique()
            adata_1hr_rem_2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_1hr_type_4_remodeling_std_mul_2_sample_num_2.h5ad')
            adata_1hr_rem_2.var_names_make_unique()
            adata_1hr_rem_3 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_1hr_type_4_remodeling_std_mul_2_sample_num_3.h5ad')
            adata_1hr_rem_3.var_names_make_unique()
            
            adata_4hr_rem_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_4hr_type_4_remodeling_std_mul_2_sample_num_1.h5ad')
            adata_4hr_rem_1.var_names_make_unique()
            adata_4hr_rem_2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_4hr_type_4_remodeling_std_mul_2_sample_num_2.h5ad')
            adata_4hr_rem_2.var_names_make_unique()
            adata_4hr_rem_3 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_4hr_type_4_remodeling_std_mul_2_sample_num_3.h5ad')
            adata_4hr_rem_3.var_names_make_unique()

            adata_D3_1_rem_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D3_1_type_4_remodeling_std_mul_2_sample_num_1.h5ad')
            adata_D3_1_rem_1.var_names_make_unique()
            adata_D3_1_rem_2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D3_1_type_4_remodeling_std_mul_2_sample_num_2.h5ad')
            adata_D3_1_rem_2.var_names_make_unique()
            adata_D3_1_rem_3 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D3_1_type_4_remodeling_std_mul_2_sample_num_3.h5ad')
            adata_D3_1_rem_3.var_names_make_unique()

            adata_D3_3_rem_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D3_3_type_4_remodeling_std_mul_2_sample_num_1.h5ad')
            adata_D3_3_rem_1.var_names_make_unique()
            adata_D3_3_rem_2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D3_3_type_4_remodeling_std_mul_2_sample_num_2.h5ad')
            adata_D3_3_rem_2.var_names_make_unique()
            adata_D3_3_rem_3 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D3_3_type_4_remodeling_std_mul_2_sample_num_3.h5ad')
            adata_D3_3_rem_3.var_names_make_unique()

            adata_D7_2_rem_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D7_2_type_4_remodeling_std_mul_2_sample_num_1.h5ad')
            adata_D7_2_rem_1.var_names_make_unique()
            adata_D7_2_rem_2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D7_2_type_4_remodeling_std_mul_2_sample_num_2.h5ad')
            adata_D7_2_rem_2.var_names_make_unique()
            adata_D7_2_rem_3 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D7_2_type_4_remodeling_std_mul_2_sample_num_3.h5ad')
            adata_D7_2_rem_3.var_names_make_unique()

            adata_D7_3_rem_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D7_3_type_4_remodeling_std_mul_2_sample_num_1.h5ad')
            adata_D7_3_rem_1.var_names_make_unique()
            adata_D7_3_rem_2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D7_3_type_4_remodeling_std_mul_2_sample_num_2.h5ad')
            adata_D7_3_rem_2.var_names_make_unique()
            adata_D7_3_rem_3 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Perturbed_adatas/adata_perturbed_D7_3_type_4_remodeling_std_mul_2_sample_num_3.h5ad')
            adata_D7_3_rem_3.var_names_make_unique()

            return {
                'Sham_1': adata_Sham_1,
                '1hr_rem_1': adata_1hr_rem_1,
                '1hr_rem_2': adata_1hr_rem_2,
                '1hr_rem_3': adata_1hr_rem_3,
                
                '4hr_rem_1': adata_4hr_rem_1,
                '4hr_rem_2': adata_4hr_rem_2,
                '4hr_rem_3': adata_4hr_rem_3,

                'D3_1_rem_1': adata_D3_1_rem_1,
                'D3_1_rem_2': adata_D3_1_rem_2,
                'D3_1_rem_3': adata_D3_1_rem_3,

                'D3_3_rem_1': adata_D3_3_rem_1,
                'D3_3_rem_2': adata_D3_3_rem_2,
                'D3_3_rem_3': adata_D3_3_rem_3,

                'D7_2_rem_1': adata_D7_2_rem_1,
                'D7_2_rem_2': adata_D7_2_rem_2,
                'D7_2_rem_3': adata_D7_2_rem_3,

                'D7_3_rem_1': adata_D7_3_rem_1,
                'D7_3_rem_2': adata_D7_3_rem_2,
                'D7_3_rem_3': adata_D7_3_rem_3,
            }
        elif dataset == "merfish_mouse_heart":
            # adata_4wk = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/merfish_mouse_heart/adata4wk.h5ad')
            # adata_4wk.var_names_make_unique()
            # adata_24wk = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/merfish_mouse_heart/adata24wk.h5ad')
            # adata_24wk.var_names_make_unique()
            # adata_90wk = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/merfish_mouse_heart/adata90wk.h5ad')
            # adata_90wk.var_names_make_unique()

            # adata_4wk.obsm['spatial'] = adata_4wk.obsm['spatial'].astype('float32')
            # adata_24wk.obsm['spatial'] = adata_24wk.obsm['spatial'].astype('float32')
            # adata_90wk.obsm['spatial'] = adata_90wk.obsm['spatial'].astype('float32')

            adata_4wk_slice_0_donor_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/merfish_mouse_heart/adata_4wk_slice_0_donor_1.h5ad')
            adata_4wk_slice_0_donor_1.var_names_make_unique()

            adata_24wk_slice_0_donor_10 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/merfish_mouse_heart/adata_24wk_slice_0_donor_10.h5ad')
            adata_24wk_slice_0_donor_10.var_names_make_unique()

            return {
                '4wk_slice_0_donor_1': adata_4wk_slice_0_donor_1,
                '24wk_slice_0_donor_10': adata_24wk_slice_0_donor_10,
                # '90wk': adata_90wk,
            }
        elif dataset == "mouse_AD_raw":
            samples = ['B02_D1', 'B02_E1', 'B03_C2', 'B03_D2', 'N02_C1', 'N02_D1', 'N03_C2', 'N03_D2', 'B07_C2', 'N07_C1', 'B06_E1', 'N06_D2', 'B04_D1', 'B04_E1', 'B05_D2', 'B05_E2', 'N04_D1', 'N04_E1', 'N05_C2', 'N05_D2']
            adatas_map = {}
            for sample in samples:
                adatas_map[sample] = sc.read(f'/home/nuwaisir/Corridor/Samee_sir_lab/Data/Mouse_AD/raw_adatas/adata_raw_{sample}.h5ad')

            return adatas_map
        elif dataset == "mouse_AD_log":
            samples = ['B02_D1', 'B02_E1', 'B03_C2', 'B03_D2', 'N02_C1', 'N02_D1', 'N03_C2', 'N03_D2', 'B07_C2', 'N07_C1', 'B06_E1', 'N06_D2', 'B04_D1', 'B04_E1', 'B05_D2', 'B05_E2', 'N04_D1', 'N04_E1', 'N05_C2', 'N05_D2']
            adatas_map = {}
            for sample in samples:
                adatas_map[sample] = sc.read(f'/home/nuwaisir/Corridor/Samee_sir_lab/Data/Mouse_AD/log_adatas/adata_log_{sample}.h5ad')

            return adatas_map
        
        elif dataset == 'mouse_AD_STARMAP_PLUS':
            adata_13months_control_rep1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/Mouse_AD_STARMAP_PLUS/adatas/adata_13months_control_replicate_1_9498.h5ad')
            adata_13months_control_rep2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/Mouse_AD_STARMAP_PLUS/adatas/adata_13months_control_replicate_2_11351.h5ad')
            adata_13months_control_rep1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/Mouse_AD_STARMAP_PLUS/adatas/adata_13months_disease_replicate_1_9494.h5ad')
            adata_13months_disease_rep2 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/Mouse_AD_STARMAP_PLUS/adatas/adata_13months_disease_replicate_2_11346.h5ad')

            return {
                'adata_13months_control_rep1': adata_13months_control_rep1,
                'adata_13months_control_rep2': adata_13months_control_rep2,
                'adata_13months_disease_rep1': adata_13months_control_rep1,
                'adata_13months_disease_rep2': adata_13months_disease_rep2,
            }
        elif dataset == 'King_gene_exp_turned_off':
            adata_Sham_1 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_Sham_1.h5ad')
            adata_d3_3 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_D3_3.h5ad')
            adata_d3_3_adata_D3_3_Prelid1_Arg1_S100a8_turned_off = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_D3_3_Prelid1_Arg1_S100a8_turned_off.h5ad')
            adata_d3_3_adata_D3_3_Got2_turned_off = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_D3_3_Got2_turned_off.h5ad')
            adata_d3_3_adata_D3_3_Got2_turned_off = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_D3_3_Got2_turned_off.h5ad')
            adata_d3_3_adata_D3_3_Pdcd1lg2_turned_off = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_D3_3_Pdcd1lg2_turned_off.h5ad')
            adata_d3_3_adata_D3_3_Pdcd1lg2_turned_off_zone_IZ = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_D3_3_Pdcd1lg2_turned_off_zone_IZ.h5ad')
            adata_d3_3_top_10_genes_off = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_D3_3_Pdcd1lg2_Entpd3_Gpc3_Card9_Il6_Olfr57_Trem1_1700113A16Rik_Pcdhb11_Sh3gl3_Pdcd1lg2_turned_off.h5ad')
            adata_d3_1_top_10_genes_off = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_D3_1_Pdcd1lg2_Entpd3_Gpc3_Card9_Il6_Olfr57_Trem1_1700113A16Rik_Pcdhb11_Sh3gl3_Pdcd1lg2_turned_off.h5ad')
            adata_Sham_1_hvg_2000_top_10_genes_off = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_Sham_1_hvg_2000_Pdcd1lg2_Entpd3_Gpc3_Card9_Il6_Olfr57_Trem1_1700113A16Rik_Pcdhb11_Sh3gl3_Pdcd1lg2_turned_off.h5ad')
            adata_Sham_1_hvg_2000 = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_Sham_1_hvg_2000.h5ad')
            adata_d3_1_hvg_2000_top_10_genes_off = sc.read('/home/nuwaisir/Corridor/Samee_sir_lab/Data/King/Fixed_adatas/adata_D3_1_hvg_2000_Pdcd1lg2_Entpd3_Gpc3_Card9_Il6_Olfr57_Trem1_1700113A16Rik_Pcdhb11_Sh3gl3_Pdcd1lg2_turned_off.h5ad')

            return {
                'adata_Sham_1': adata_Sham_1,
                'adata_Sham_1_hvg_2000': adata_Sham_1_hvg_2000,
                'adata_d3_3': adata_d3_3,
                'adata_d3_3_Prelid1_Arg1_S100a8_turned_off': adata_d3_3_adata_D3_3_Prelid1_Arg1_S100a8_turned_off,
                'adata_d3_3_Got2_turned_off': adata_d3_3_adata_D3_3_Got2_turned_off,
                'adata_d3_3_Pdcd1lg2_turned_off': adata_d3_3_adata_D3_3_Pdcd1lg2_turned_off,
                'adata_d3_3_Pdcd1lg2_turned_off_zone_IZ': adata_d3_3_adata_D3_3_Pdcd1lg2_turned_off_zone_IZ,
                'adata_Sham_1_hvg_2000_top_10_genes_off': adata_Sham_1_hvg_2000_top_10_genes_off,
                'adata_d3_1_hvg_2000_top_10_genes_off': adata_d3_1_hvg_2000_top_10_genes_off,
                'adata_d3_3_top_10_genes_off': adata_d3_3_top_10_genes_off,
                'adata_d3_1_top_10_genes_off': adata_d3_1_top_10_genes_off,
            }

        
        else:
            print("Can't find the specified dataset!")
            return []
