import numpy as np
import os
import shutil
import matplotlib.pyplot as plt
from runFEA import runFEA
from getIMG import getIMG
from CalSF import CalSF
from Calculation import Calculation
from genFEA import genFEA
from genTar import genTar
import time

plt.style.use('seaborn-whitegrid')

fname = os.getcwd() + "/"

if os.path.isdir("05_mesh_00") == False:
    shutil.copytree("05_mesh_ref", "05_mesh_00")

m_a = 1585 # 1585
m_b = 7.41
m_af = 10262
m_bf = 28.882 # 28.882
m_as = 271
m_bs = 4.08
m_afs = 2160
m_bfs = 11.25
m_eta = 400 # 400
m_E = 1.0e5

mat_para_opt = np.array([m_a, m_b, m_af, m_bf, m_as, m_bs, m_afs, m_bfs])
mat_para_fix = np.array([m_eta, m_E])

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
    if r_disp_max < 0.01:
        break
    else:
        # f_ind = 3
        print("++++++++++++++++++++++++++++++++++++++")
        print("The " + str(f_ind) + " generation: ")

        runFEA(mat_para_opt, f_ind, mat_para_fix, finalflag)
        vals = genFEA(fname, f_ind, np.array([1]), p_info)
        imageData = getIMG(imagePath, file_index, vals, p_info)
        cal_vals = CalSF(fname, f_ind, vals, imageData)
        genTar(fname, f_ind, vals, cal_vals)
        r_disp_max = cal_vals['r_disp_max']

print("++++++++++++++++++++++++++++++++++++++")
print("The final " + str(f_ind) + " generation: ")
finalflag = True
runFEA(mat_para_opt, f_ind, mat_para_fix, finalflag)
vals = genFEA(fname, f_ind, file_sub, p_info)
cal_vals = Calculation(fname, f_ind, vals, imageData)

k = '%02d' % (f_ind + 1)
shutil.rmtree("05_mesh_" + k, ignore_errors=True)
print("Elapsed: " + str(time.time() - t))

print()



