import numpy as np
import os
import shutil
from getIMG import getIMG
from genFEA import genFEA

fname = os.getcwd() + "/"

if os.path.isdir("05_mesh_00") == False:
    shutil.copytree("05_mesh_ref", "05_mesh_00")

lv_p_gindex = np.array([3769,2138,1591])
rv_p_gindex = np.array([7289,7421,1548])
epi_p_gindex = np.array([6314,4259,1353])
p_info = [lv_p_gindex, rv_p_gindex, epi_p_gindex]

fname = os.getcwd() + "/"
imagePath = "morph/"
file_index = ["01_RR_70", "02_RR_80", "03_RR_90", "04_RR_99"]
# file_sub = np.array([0.25,0.35,0.75,1])
file_sub = np.linspace(0.05, 1.0, num=20, endpoint=True)

f_ind = 2

vals = genFEA(fname, f_ind, file_sub, p_info)

imageData = getIMG(imagePath, file_index, vals, p_info)

lv_vlm_img = imageData['lv_vlm_img']
rv_vlm_img = imageData['rv_vlm_img']
lv_vlm_cur = vals['lv_vlm_cur']
rv_vlm_cur = vals['rv_vlm_cur']
print("LV Volume FEA: ")
print(vals['lv_vlm_dat'], vals['lv_vlm_ref'], lv_vlm_cur)
print("RV Volume FEA: ")
print(vals['rv_vlm_dat'], vals['rv_vlm_ref'], rv_vlm_cur)
print("LV Volume IMAGE: ")
print(lv_vlm_img)
print("RV Volume IMAGE: ")
print(rv_vlm_img)


print()



