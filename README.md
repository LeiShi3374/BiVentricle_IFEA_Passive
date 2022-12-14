# BiVentricle_IFEA_Passive

This repository contains the source code for our paper:

Lei Shi, Hiroo Takayama, Vijay Vedula. An Optimization Scheme to Quantify Late-diastolic Bi-ventricular Mechanics in Patients with Obstructive Hypertrophic Cardiomyopathy (HOCM) 

## Dependency

### Python packages
VTK

numpy

### svFSI
https://github.com/LeiShi3374/svFSI 

## Usage 

Use Genetic Algorithm to search for the passive material parameters for the bi-ventricular model and find the stress-free configuration at the same time. 

### Hyperparameter and variable boundary 

lb and ub defines the lower and upper boundaries of the material parameters. 

The GA algorithm accepts several hyperparameters including the maximum number of iterations, the mutation rate, the number of population (individuals), and etc. 

### Hyperelastic model 

The example used here is using the HO model (Holzapfel, Gerhard A., and Ray W. Ogden. "Constitutive modelling of arteries." Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences 466.2118 (2010): 1551-1597.) 

If other model is used, the corresponding sections in the main function and the runFEA.py function need to be revised. 

In the current example, 6 HO material parameters a, b, a_f, b_f, a_s, b_s are set for free to optimize. 

### Mesh data 

The mesh data including the fiber orientations are included in the folder "05_mesh_ref". The fiber orientations are generated using a rule-based method following the method in the paper Bayer, Jason D., et al. "A novel rule-based algorithm for assigning myocardial fiber orientation to computational heart models." Annals of biomedical engineering 40.10 (2012): 2243-2254. 

### Experimental data 

The folder "morph" contains the registered image data obtained using the morphing algorithm, serving as the real data to be captured using the model. In the current
example, the image data of 4 time points (70%, 80%, 90%, 90% RR-interval) are used which cover the deformation from the diastasis to the end-diastolic states. 

### Boundary conditions 

The pressures exerted on the LV and RV endocardium surfaces are from literature (pressure_lv_p.dat and pressure_rv_p.dat). 

The Robin boundaries are used to simulate the pericardium effects. (Pfaller, Martin R., et al. "The importance of the pericardium for cardiac biomechanics: from physiology to computational modeling." Biomechanics and modeling in mechanobiology 18.2 (2019): 503-529.)

### Objective function

The objective function to minimize is the difference of the LV and RV volume changes, and the displacements of several landmarks between the image and FEA data.   

9 landmarks are selected from the left and right endocardium surfaces and the epicardium surface. 

Currently, only Genetic Algorithm is supported in this framework, and it will be updated to also support the Bayesian optimization and MCMC. 








