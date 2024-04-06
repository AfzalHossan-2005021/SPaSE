import numpy as np
import pandas as pd
# import paste as pst
from anndata import AnnData
from typing import List, Tuple, Optional
import torch
import ot
import scipy
import matplotlib.pyplot as plt
import warnings
from .DataLoader import DataLoader
from .Preprocessor import Preprocessor
from .utils import trim_barcode, make_format, norm_and_center_coordinates, rotate, mirror, euc_dist, plot_slice_pairwise_alignment_modified, scale_coords, QC, to_dense_array, extract_data_matrix, intersect, kl_divergence_backend, jensenshannon_divergence_backend
import os
from .paste_helper_functions import match_spots_using_spatial_heuristic
import json
import scanpy as sc


class PairwiseAlign():
    def __init__(self, config):
        self.config = config
        self.dataset = config['dataset']
        self.sample_left = config['sample_left']
        self.sample_right = config['sample_right']
        self.alpha = config['alpha']
        self.dissimilarity = config['dissimilarity']
        self.init_map_scheme = config['init_map_scheme']
        self.numIterMaxEmd = config['numIterMaxEmd']
        self.numInnerIterMax = config['numInnerIterMax']
        self.use_gpu = config['use_gpu']
        self.results_path = config['results_path']
        self.config_file_name = os.path.basename(config['config_path'])
        self.sinkhorn = config['sinkhorn']
        self.lambda_sinkhorn = config['lambda_sinkhorn']
        self.cost_mat_path = f'{self.results_path}/../local_data/{self.dataset}/{self.sample_left}/cost_mat_{self.sample_left}_{self.sample_right}_{self.dissimilarity}.npy'
        os.makedirs(os.path.dirname(self.cost_mat_path), exist_ok=True)

        if config['adata_left_path'] != 'None':
            self.adata_left = sc.read(config['adata_left_path'])
            self.adata_right = sc.read(config['adata_right_path'])
        else:
            data_loader = DataLoader(config)

            dataset_map = data_loader.read_data(self.dataset)

            self.pi_low_entropy_path = f'{self.results_path}/{self.dataset}/config_{self.dataset}_{self.sample_left}_vs_{self.sample_right}_js.json/Pis/{self.dataset}_uniform_js.npy'

            self.adata_left = dataset_map[self.sample_left]
            self.adata_right = dataset_map[self.sample_right]

        scale_coords(self.adata_left, key_name='spatial')
        scale_coords(self.adata_right, key_name='spatial')

        if config['QC']:
            QC(self.adata_left)
            QC(self.adata_right)

    def pairwise_align_sinkhorn(self):
        if not torch.cuda.is_available():
            if self.use_gpu == True:
                print("Setting use_gpu to False")
            self.use_gpu = False

        if self.use_gpu:
            backend = ot.backend.TorchBackend()
        else:
            backend = ot.backend.NumpyBackend()
        pi_init = None

        if self.init_map_scheme == 'uniform':
            pi_init = None
            
        print("Calculating pi using gcg")

        # print("......... Attention!!! Need to pass in cost_mat_path! ........")
        self.config['cost_mat_path'] = self.cost_mat_path
        os.makedirs(os.path.dirname(self.cost_mat_path), exist_ok=True)

        if self.lambda_sinkhorn == 'inf':
            pi = np.ones((self.adata_left.n_obs, self.adata_right.n_obs)) / (self.adata_left.n_obs * self.adata_right.n_obs)
            if self.use_gpu:
                pi = torch.from_numpy(pi)
            return pi, -1000
        
        pi, fgw_dist = self.pairwise_align_modified(self.adata_left,self.adata_right,alpha=self.alpha,sinkhorn=self.sinkhorn,lambda_sinkhorn=self.lambda_sinkhorn,dissimilarity=self.dissimilarity,G_init=pi_init, numItermax=10000,cost_mat_path=self.cost_mat_path,return_obj=True,norm=True,verbose=False,backend=backend,use_gpu=self.use_gpu,numInnerItermax=self.numInnerIterMax)
        return pi, fgw_dist


    def pairwise_align_modified(
        self,
        sliceA: AnnData, 
        sliceB: AnnData, 
        alpha: float = 0.1, 
        dissimilarity: str = 'js', 
        sinkhorn: bool = False,
        use_rep: Optional[str] = None,
        lambda_sinkhorn: float = 1, 
        G_init = None, 
        a_distribution = None, 
        b_distribution = None, 
        norm: bool = True, 
        numItermax: int = 10000, 
        backend = ot.backend.NumpyBackend(), 
        use_gpu: bool = False, 
        return_obj: bool = False, 
        verbose: bool = False, 
        gpu_verbose: bool = True,
        cost_mat_path: Optional[str] = None,
        **kwargs) -> Tuple[np.ndarray, Optional[int]]:
        """
        Calculates and returns optimal alignment of two slices. 
        
        Args:
            sliceA: Slice A to align.
            sliceB: Slice B to align.
            alpha:  Alignment tuning parameter. Note: 0 <= alpha <= 1.
            dissimilarity: Expression dissimilarity measure: ``'kl'`` or ``'euclidean'`` or ``'jensenshannon'``.
            use_rep: If ``None``, uses ``slice.X`` to calculate dissimilarity between spots, otherwise uses the representation given by ``slice.obsm[use_rep]``.
            G_init (array-like, optional): Initial mapping to be used in FGW-OT, otherwise default is uniform mapping.
            a_distribution (array-like, optional): Distribution of sliceA spots, otherwise default is uniform.
            b_distribution (array-like, optional): Distribution of sliceB spots, otherwise default is uniform.
            numItermax: Max number of iterations during FGW-OT.
            norm: If ``True``, scales spatial distances such that neighboring spots are at distance 1. Otherwise, spatial distances remain unchanged.
            backend: Type of backend to run calculations. For list of backends available on system: ``ot.backend.get_backend_list()``.
            use_gpu: If ``True``, use gpu. Otherwise, use cpu. Currently we only have gpu support for Pytorch.
            return_obj: If ``True``, additionally returns objective function output of FGW-OT.
            verbose: If ``True``, FGW-OT is verbose.
            gpu_verbose: If ``True``, print whether gpu is being used to user.
    
        Returns:
            - Alignment of spots.

            If ``return_obj = True``, additionally returns:
            
            - Objective function output of FGW-OT.
        """
        
        # Determine if gpu or cpu is being used
        if use_gpu:
            try:
                import torch
            except:
                print("We currently only have gpu support for Pytorch. Please install torch.")
                    
            if isinstance(backend,ot.backend.TorchBackend):
                if torch.cuda.is_available():
                    if gpu_verbose:
                        print("gpu is available, using gpu.")
                else:
                    if gpu_verbose:
                        print("gpu is not available, resorting to torch cpu.")
                    use_gpu = False
            else:
                print("We currently only have gpu support for Pytorch, please set backend = ot.backend.TorchBackend(). Reverting to selected backend cpu.")
                use_gpu = False
        else:
            if gpu_verbose:
                print("Using selected backend cpu. If you want to use gpu, set use_gpu = True.")
                
        # subset for common genes
        common_genes = intersect(sliceA.var.index, sliceB.var.index)
        sliceA = sliceA[:, common_genes]
        sliceB = sliceB[:, common_genes]

        # Backend
        nx = backend    
        
        # Calculate spatial distances
        coordinatesA = sliceA.obsm['spatial'].copy()
        coordinatesA = nx.from_numpy(coordinatesA)
        coordinatesB = sliceB.obsm['spatial'].copy()
        coordinatesB = nx.from_numpy(coordinatesB)
        
        if isinstance(nx,ot.backend.TorchBackend):
            coordinatesA = coordinatesA.float()
            coordinatesB = coordinatesB.float()

        D_A = ot.dist(coordinatesA,coordinatesA, metric='euclidean')
        D_B = ot.dist(coordinatesB,coordinatesB, metric='euclidean')

        if isinstance(nx,ot.backend.TorchBackend) and use_gpu:
            D_A = D_A.cuda()
            D_B = D_B.cuda()
        
        # Calculate expression dissimilarity
        A_X, B_X = nx.from_numpy(to_dense_array(extract_data_matrix(sliceA,use_rep))), nx.from_numpy(to_dense_array(extract_data_matrix(sliceB,use_rep)))

        if isinstance(nx,ot.backend.TorchBackend) and use_gpu:
            A_X = A_X.cuda()
            B_X = B_X.cuda()

        if os.path.exists(cost_mat_path):
            print("Loading cost matrix from file system...")
            M = np.load(cost_mat_path)
        else:
            print("cost_mat_path does not exist.")
            if dissimilarity.lower()=='euclidean' or dissimilarity.lower()=='euc':
                M = ot.dist(A_X,B_X)
            elif dissimilarity.lower()=='kl':
                s_A = A_X + 0.01
                s_B = B_X + 0.01
                M = kl_divergence_backend(s_A, s_B)
            elif dissimilarity.lower()=='js' or dissimilarity.lower()=='jensenshannon':
                s_A = A_X + 0.01
                s_B = B_X + 0.01
                M = jensenshannon_divergence_backend(s_A, s_B)
            np.save(cost_mat_path, M)
        M = nx.from_numpy(M)

        if isinstance(nx,ot.backend.TorchBackend) and use_gpu:
            M = M.cuda()
        
        # init distributions 
        if a_distribution is None:
            a = nx.ones((sliceA.shape[0],))/sliceA.shape[0]
        else:
            a = nx.from_numpy(a_distribution)
            
        if b_distribution is None:
            b = nx.ones((sliceB.shape[0],))/sliceB.shape[0]
        else:
            b = nx.from_numpy(b_distribution)

        if isinstance(nx,ot.backend.TorchBackend) and use_gpu:
            a = a.cuda()
            b = b.cuda()
        
        if norm:
            D_A /= nx.min(D_A[D_A>0])
            D_B /= nx.min(D_B[D_B>0])
        
        # Run OT
        if G_init is not None:
            G_init = nx.from_numpy(G_init)
            if isinstance(nx,ot.backend.TorchBackend):
                G_init = G_init.float()
                if use_gpu:
                    G_init.cuda()
        
        assert(sinkhorn == 1)

        pi, logw = self.my_fused_gromov_wasserstein_gcg(M, D_A, D_B, a, b, lambda_sinkhorn=lambda_sinkhorn, G_init = G_init, loss_fun='square_loss', alpha= alpha, log=True, numItermax=numItermax,verbose=verbose, use_gpu = use_gpu, **kwargs)
        
        pi = nx.to_numpy(pi)
        obj = nx.to_numpy(logw['fgw_dist'])
        if isinstance(backend,ot.backend.TorchBackend) and use_gpu:
            torch.cuda.empty_cache()

        if return_obj:
            return pi, obj
        return pi
    
    
    def my_fused_gromov_wasserstein_gcg(self, M, C1, C2, p, q, lambda_sinkhorn=1, G_init = None, loss_fun='square_loss', alpha=0.5, armijo=False, log=False,numItermax=200, use_gpu = False, **kwargs):
        """
        Adapted fused_gromov_wasserstein with the added capability of defining a G_init (inital mapping).
        Also added capability of utilizing different POT backends to speed up computation.
        
        For more info, see: https://pythonot.github.io/gen_modules/ot.gromov.html
        """
        # print("Inside my_fused_gromov_wasserstein_gcg")
        # print(f'alpha: {alpha}')

        p, q = ot.utils.list_to_array(p, q)

        p0, q0, C10, C20, M0 = p, q, C1, C2, M
        # print('C10')
        # print(C10)
        nx = ot.backend.get_backend(p0, q0, C10, C20, M0)

        constC, hC1, hC2 = ot.gromov.init_matrix(C1, C2, p, q, loss_fun)

        if G_init is None:
            G0 = p[:, None] * q[None, :]
        else:
            G0 = (1/nx.sum(G_init)) * G_init
            if use_gpu:
                G0 = G0.cuda()

        # print('hC1:')
        # print(hC1)

        # print('hC2:')
        # print(hC2)

        def f(G):
            return ot.gromov.gwloss(constC, hC1, hC2, G)

        def df(G):
            return ot.gromov.gwggrad(constC, hC1, hC2, G)

        if log:
            # print('doing gcg')
            # print((1 - alpha) * M)
            res, log = ot.optim.gcg(p, q, M, lambda_sinkhorn, alpha, f, df, G0, log=True, **kwargs)
            # res, log = ot.gromov.cg(p, q, (1 - alpha) * M, alpha, f, df, G0, armijo=armijo, C1=C1, C2=C2, constC=constC, log=True, **kwargs)

            fgw_dist = log['loss'][-1]

            log['fgw_dist'] = fgw_dist
            # log['u'] = log['u']
            # log['v'] = log['v']
            return res, log

        else:
            # return ot.gromov.cg(p, q, (1 - alpha) * M, alpha, f, df, G0, armijo=armijo, C1=C1, C2=C2, constC=constC, **kwargs)
            # print('pi before gcg')
            pi = ot.optim.gcg(p, q, M, lambda_sinkhorn, alpha, f, df, G0, log=False, **kwargs)
            # print('pi after gcg')
            return pi, -1