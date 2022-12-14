'''This code is designed to find both the stress-free configuration
and the optimized material parameters for the passive process
from the diastasis to the end-diastolic states
HO model is used as the constitutive model and GA is used for
the IFEA process
by Lei Shi at Vedula Lab in Columbia University 2022
'''

import numpy as np
import os
import shutil
from runFEA import runFEA
from getIMG import getIMG
from CalSF import CalSF
from Calculation import Calculation
from genFEA import genFEA
from genTar import genTar
from scipy.optimize import least_squares
import GeAlgo as ga
import time

def Simulations(x):

    if os.path.isdir("05_mesh_00") == False:
        shutil.copytree("05_mesh_ref", "05_mesh_00")

    # Get material parameters
    m_a = 10 ** x[0]
    m_b = x[1]
    m_af = 10 ** x[2]
    m_bf = x[3]
    m_as = 10 ** x[4]
    m_bs = x[5]
    m_afs = 2160
    m_bfs = 11.25
    m_eta = 400 # 400
    m_E = 1.0e5

    mat_para_opt = np.array([m_a, m_b, m_af, m_bf, m_as, m_bs, m_afs, m_bfs])
    mat_para_fix = np.array([m_eta, m_E])

    # Node index for selected landmarks
    lv_p_gindex = np.array([3769,2138,1591])
    rv_p_gindex = np.array([7289,7421,1548])
    epi_p_gindex = np.array([6314,4259,1353])
    p_info = [lv_p_gindex, rv_p_gindex, epi_p_gindex]

    fname = os.getcwd() + "/"
    imagePath = "morph/"
    file_index = ["01_RR_70", "02_RR_80", "03_RR_90", "04_RR_99"]
    file_sub = np.array([0.25,0.5,0.75,1])

    r_disp_max = 1
    finalflag = False

    t = time.time()

    for f_ind in range(3):
        # f_ind = 3
        print("++++++++++++++++++++++++++++++++++++++")
        print("The " + str(f_ind) + " generation: ")

        runFEA(mat_para_opt, f_ind, mat_para_fix, finalflag)
        vals = genFEA(fname, f_ind, np.array([1]), p_info)
        imageData = getIMG(imagePath, file_index, vals, p_info)
        cal_vals = CalSF(fname, f_ind, vals, imageData)
        genTar(fname, f_ind, vals, cal_vals)
        r_disp_max = cal_vals['r_disp_max']
        if r_disp_max < 0.01:
            break


    print("++++++++++++++++++++++++++++++++++++++")
    print("The final " + str(f_ind) + " generation: ")
    finalflag = True
    runFEA(mat_para_opt, f_ind, mat_para_fix, finalflag)
    vals = genFEA(fname, f_ind, file_sub, p_info)
    cal_vals = Calculation(fname, f_ind, vals, imageData)

    k = '%02d' % (f_ind + 1)
    shutil.rmtree("05_mesh_" + k, ignore_errors=True)
    print("Elapsed: " + str(time.time() - t))

    fit_err = cal_vals['fit_err']
    if np.isnan(fit_err):
        fit_err = 100

    return fit_err

tinit = time.time()

# GA #################################################
# lower and upper boundaries for parameters
# a, b, a_f, b_f, a_s, b_s NOTE: a, a_f, a_s are log10-ed to
# avoid order effects
lb = np.array([2,1,3,1,3,1])
ub = np.array([5,40,6,40,6,40])

algorithm_parameter = {'max_num_iteration': 30,
                   'population_size':4,
                   'mutation_probability':0.1,
                    'mutation_change_generation': 5,
                    'mutation_change_factor': 1.3}

r = ga.GeAlgo(Simulations, lb, ub, algorithm_parameter)
# GA #################################################
# Bayesian ###########################################

# Bayesian ###########################################

print("Total Elapsed: " + str(time.time() - tinit))

print()



