<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Don't forget to give the project a star!
*** Thanks again! Now go create something AMAZING! :D
-->



<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
<!-- [![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url] -->



<!-- PROJECT LOGO -->
<!-- <br />
<div align="center">
  <a href="https://github.com/othneildrew/Best-README-Template">
    <img src="images/logo.png" alt="Logo" width="80" height="80">
  </a>

  <h3 align="center">Best-README-Template</h3>

  <p align="center">
    An awesome README template to jumpstart your projects!
    <br />
    <a href="https://github.com/othneildrew/Best-README-Template"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/othneildrew/Best-README-Template">View Demo</a>
    ·
    <a href="https://github.com/othneildrew/Best-README-Template/issues/new?labels=bug&template=bug-report---.md">Report Bug</a>
    ·
    <a href="https://github.com/othneildrew/Best-README-Template/issues/new?labels=enhancement&template=feature-request---.md">Request Feature</a>
  </p>
</div> -->




<!-- ## Table of Contents

<ol>
<li>
    <a href="#SPaSE">SPaSE</a>
    <ul>
    <li><a href="#built-with">Built With</a></li>
    </ul>
</li>
<li>
    <a href="#getting-started">Getting Started</a>
    <ul>
    <li><a href="#prerequisites">Prerequisites</a></li>
    <li><a href="#installation">Installation</a></li>
    </ul>
</li>
<li><a href="#usage">Usage</a></li>
<li><a href="#roadmap">Roadmap</a></li>
<li><a href="#contributing">Contributing</a></li>
<li><a href="#license">License</a></li>
<li><a href="#contact">Contact</a></li>
<li><a href="#acknowledgments">Acknowledgments</a></li>
</ol> -->



<!-- ABOUT THE PROJECT -->
## SPaSE
SPaSE (Spatially-resolved Pathology ScorE), designed to quantify pathological effects within an ST tissue sample incorporating an optimal transport problem formulation between the pathologically impacted and control reference ST samples, considering both gene expression and spatial spot locations.


<p align="right">(<a href="#readme-top">back to top</a>)</p>



### Built With

[![Python][python-img]][python-url]


<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- GETTING STARTED -->
## Getting Started

### Installation


1. Clone the repo
   ```
   git clone https://github.com/Nuwaisir-1998/SPaSE.git
   ```
2. Install the spase conda environment:
   ```
   conda env create -f environment.yml
   ```

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Prerequisite for Running the Method
We utilized the POT (Python OT) package to calculate the optimal transport plan, specifically employing the 'sinkhorn_log' method to address the entropic regularization optimal transport problem. To implement this change, navigate to the directory where your conda environments are stored. Proceed to the directory containing all the Python packages for the 'spase' conda environment. Within that directory, locate the 'ot' directory, where you'll find the 'bregman.py' file. Adjust the initialization of the method variable to 'sinkhorn_log' in this file. In my case, the path for 'bregman.py' is ```/home/nuwaisir/miniconda3/envs/paste/lib/python3.8/site-packages/ot/bregman.py```. The modified function definition should resemble the following:
```
def sinkhorn(a, b, M, reg, method='sinkhorn_log', numItermax=1000,
             stopThr=1e-9, verbose=False, log=False, warn=True,
             **kwargs):
```

## Usage

1. The anndata objects for both healthy and diseased samples are accessible via the following link: [Data link](https://drive.google.com/drive/folders/1gXsxIYXrSixPAb6l09L4IYycLqYLLngU?usp=sharing). Below are descriptions of the files:
    - ```adata_Sham_1.h5ad```: healthy sample
    - ```adata_1hr.h5ad```: 1-hour post-MI sample
    - ```adata_4hr.h5ad```: 4-hour post-MI sample
    - ```adata_D3_1.h5ad```: 3-day post-MI sample (replicate 1)
    - ```adata_D3_3.h5ad```: 3-day post-MI sample (replicate 3)
    - ```adata_D7_2.h5ad```: 7-day post-MI sample (replicate 2)
    - ```adata_D7_3.h5ad```: 7-day post-MI sample (replicate 3)

2. Place all the files inside ```Data/King/Fixed_adatas/```
3. Naviage to ```Workspace/SPaSE/src/```
4. The following command will reproduce the optimal transport plans and the pathological scores for each healthy disease sample pair. ```python rep_run.py```
5. If you want to run SPaSE on two samples of your choice, then inside ```run.py```, set the values of the following variable:
    - dataset [type: string]
    - adata_left_path [type: string; the path of the .h5ad file that would be used as reference]
    - adata_right_path [type: string; the path of the .h5ad file that would be used as target]
    - adata_to_be_synthesized_path = [type: string; the path of the .h5ad file representing reference sample that would be decomposed for synthesization, you can assign this to adata_left_path if that is the sample you want to synthesize for calculating the null distribution]
    - sample_left [type: string; name of the reference sample]
    - sample_right [type: string; name of the target sample]
    - alpha [type: double; hyperparameter]
    - lambda_sinkhorn [type: double; hyperparameter]
6. The results will be stored inside ```Workspace/SPaSE/results/``` under the name you provide in the ```dataset``` variable.


<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



[python-img]: https://www.python.org/static/img/python-logo.png
[python-url]: https://www.python.org/
