# Statistical Inference With Betti-0 Curves

This repository provides examples of topology-inspired statistical inference using Betti-0 curves proposed in the article [Revisiting Abnormalities in Brain Network Architecture Underlying Autism Using Topology-Inspired Statistical Inference](https://doi.org/10.1089/brain.2018.0604).

The data directory contains gray matter density data for 49 subjects with autism spectrum disorders (ASD) and 49 typically developing control (TDC) subjects. The data is divided into three csv files: `7266_Raw_Densities.csv`, `ICN_ROI_Key.csv`, and `ICN_SUB_Densities.csv`.

`7266_Raw_Densities.csv`: Contains the gray matter densities from 7266 regions of interest uniformly distributed through the gray matter in the whole brain.

`ICN_ROI_Key.csv`: Contains meta-data of the ROIs underlying the salience network (SN), the executive control network (ECN), and the default mode network (DMN). These ROIs are obtained using scMRI analysis described in the 2010 article [scMRI Reveals Large-Scale Brain Network Abnormalities in Autism](https://doi.org/10.1371/journal.pone.0049172) by Zielinski et al.

`ICN_SUB_Densities.csv` Contains the gray matter densities for ROIs underlying SN, ECN, and DMN for ASD and TDC subjects.


There are three notebooks: `Examples_Betti_Curve_Plots.ipynb`, `Examples_Permutation_Test.ipynb`, and `Examples_Bootstrap_Test.ipynb`, along with a python script file `helper_functions.py`.

`helper_functions.py`: Provides functions to load the data, construct structural correlation graphs, compute minimum spanning trees and obtain Betti-0 curves.

`Examples_Betti_Curve_Plots.ipynb`: Describes how to load the data, construct structural correlation graphs (Global-SCG, SN-SCG, ECN-SCG, and DMN-SCG), and plot the Betti-0 curves.

`Examples_Permutation_Test.ipynb`: Describes the steps of permutation test on the ASD and TDC groups using Betti-0 curves.

`Examples_Bootstrap_Test.ipynb`: Describes the steps of bootstrap test on the ASD and TDC groups using Betti-0 curves.

<hr>

### Requirements
The code is written for Python 3.6 and has been tested with the following package configuration:

```
numpy==1.18.1
scipy==1.4.1              # The scipy.sparse.csgraph module is used in the helper_functions.py
matplotlib==3.1.3
```

The code should work just the same for python 2.x and/or with older versions of these packages (although I haven't tested it).

<hr>

If you find this code useful, please cite our Brain Connectivity article:

```
@Article{PalandeJoseZielinski2019,
	author   = {Palande, Sourabh and Jose, Vipin and Zielinski, Brandon and Anderson, Jeffrey and Fletcher, P. Thomas and Wang, Bei},
	journal  = {Brain Connectivity},
	title    = {Revisiting abnormalities in brain network architecture underlying autism using topology-inspired statistical inference.},
	year     = {2019},
	number   = {1},
	pages    = {13--21},
	volume   = {9}
}
```
