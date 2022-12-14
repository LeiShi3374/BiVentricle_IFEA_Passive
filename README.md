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

lb and ub defines the lower and upper boundaries of the material parameters. 

The example used here is using the HO model (Holzapfel, Gerhard A., and Ray W. Ogden. "Constitutive modelling of arteries." Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences 466.2118 (2010): 1551-1597.) 

6 material parameters a, b, a_f, b_f, a_s, b_s are set for free to optimize. 







