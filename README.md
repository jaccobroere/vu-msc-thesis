# MSc Thesis: Impulse response analysis using spatio-temporal models
*Author: Jacco Broere*

This repository contains all the necessary code for replicating the results in my MSc thesis.

## Using Kubernetes to determine the best lambda for a specific design
The following steps are taken to determine the best lambda for a specific design:
1. Create a folder called k8s_export in the root of the project, if it does not exist yet.
2. Adjust the simulation parameters in the file scripts/simulation_detlam.sh
3. Add absolute path to your /out folder in the determine_lambda.yml file
4. Rebuild the jaccusaurelius/vu-msc-thesis:kube image from the project root directory using the following command:
```docker build -t jaccusaurelius/vu-msc-thesis:kube .```
5. Run the scripts/k8s/run_kubernetes_detlam.sh script

## Using Kubernetes to run the simulation experiments with the best lambdas


