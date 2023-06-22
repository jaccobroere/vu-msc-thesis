# MSc Thesis: Graph-guided Fused Lasso approaches to the estimation of spatio-temporal autoregressions.
*Author: Jacco Broere*

This repository contains all the necessary code for replicating the results in my MSc thesis.

## 1. Structure of the repository
The main structure of the repository is as follows:

### 1.1 `src/`
The `src/` folder contains all source code that is used for computation and plotting. It consists of four subfolders:
1. `src/compute` contains all the source code related to estimation of the models and the computation of necessary preprocessing, hyperparameter optimzation etc.
2. `src/plotting` contains all the source code related to plotting the results and generating the visuals used in the thesis
3. `src/simulation` contains all the source code related to simulating the data for the simulation study in the thesis (Note that Design 1, 2, and 3 in the thesis text correspond to B, C, and D in this repository, respectively)
4. `src/playgrounds` this contains no essential code, and was mostly used to quickly explore certain ideas and test things

Note that in many of the scripts, the code makes use of an environment variable `PROJ_DIR` that is set to the root of the project. This is done to make the code portable. You should set this environment variable to run the code as is. Or you should replace the lines that make use of this with an explicit path to the root of the project. Note that you can also use the Docker image that I provide later on in this README, which already has this environment variable set.

Moreover, Note that many of the scripts import some sort of utils function at the top of the file. These utils functions are often located in some child folder or in the same folder directly.

### 1.2 `out/`
This folder is used to store the output of the estimated models, tables, figures, plots etc. Essentially, it covers all output of models.

### 1.3 `data/`
The `data/` folder is meant to store the data related to simulation, and it will also contain the precalculated matrices, such as the the penalty matrix $\tilde{\bm{D}}$, its inverse and the graphs pertaining to the different simulation designs.

### 1.4 `scripts/`
The `scripts/` folder contains all the bash scripts, and Kubernetes (k8s) configuration files that can be used to run the models and pipelines in succession for multiple designs in succession. The k8s files can be used to parallellize the computation of the simulation study over multiple cores or machines. This especially beneficial on strong machines or in the cloud (AWS, GCP, Azure etc.).

### 1.5 `docker/`
This folder contains the Dockerfile that is used to build the Docker image that is used to run the models and pipelines, such that it can be run consistently and in the same way on every machine. This is especially handy to manage package and software dependencies across machines. Note that the main Docker image to run the models is available on Dockerhub (https://hub.docker.com/repository/docker/jaccusaurelius/vu-msc-thesis/general). The `jaccusaurelius/vu-msc-thesis:workspace` is the image that was used for most of the computations in the thesis. The `Dockerfile` in this folder can be used to reconstruct a newer version of this image.

### 1.6 `admm_src_zhu/`
This is a copy of the source code that was written by Zhu et al. (2017), we use their Augmented ADMM algorithm to solve the generalized LASSO problems as part of the novel GF-SPLASH estimators we developed.

```
Zhu, Y. (2017). An Augmented ADMM Algorithm With Application to the Generalized Lasso Problem. Journal of Computational and Graphical Statistics, 26(1), 195–204. https://doi.org/10.1080/10618600.2015.1114491
```

### 1.7 `archive/`
This folder is used to store old and unused code that was written during development and testing of the code for this thesis. It is mostly very messy, however, there it can be nice to keep as there is some valuable code in there that might be useful when reviving certain interesting ideas in the future, which I wasn't able to completely or succesfully implement in time for the thesis.



## Using Kubernetes to determine the best lambda for a specific design
The following steps are taken to determine the best lambda for a specific design:
1. Create a folder called k8s_export in the root of the project, if it does not exist yet.
2. Adjust the simulation parameters in the file scripts/simulation_detlam.sh
3. Add absolute path to your /out folder in the determine_lambda.yml file
4. Rebuild the jaccusaurelius/vu-msc-thesis:kube image from the project root directory using the following command:
```docker build -t jaccusaurelius/vu-msc-thesis:kube .```
5. Run the scripts/k8s/run_kubernetes_detlam.sh script

## Using Kubernetes to run the simulation experiments with the best lambdas


