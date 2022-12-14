## This function provides a method to do optimization of the material parameters
## based on the current configuration and image configuration

import numpy as np

def Calculation(fname, f_ind, vals, imageData):
    ####################### Calculation  #######################################################

    point0d = vals['point0d']
    # point0 = vals['point0']
    point1R = vals['point1R']

    if f_ind > 0:
        k = '%02d' % (f_ind - 1)
        fname_result_pre = "05_mesh_" + k + "1/"
    k = '%02d' % (f_ind)
    fname_result = "05_mesh_" + k + "result/"

    r_disp = np.zeros(len(point0d))
    r_pos = point1R - point0d
    for i in range(len(r_pos)):
        r_disp[i] = np.linalg.norm(r_pos[i, :])
    r_disp_max = np.max(r_disp)
    np.save(fname + fname_result + "r_disp.npy", r_pos)
    print("Max disp diff: " + str(r_disp_max))

    # beta = 1
    # if "fname_result_pre" in globals():
    #     r_pos_1 = np.load(fname + fname_result_pre + "r_disp.npy")
    #     beta = -beta * np.tensordot(r_pos_1, (r_pos - r_pos_1), axes=2) / np.tensordot((r_pos - r_pos_1),
    #                                                                                    (r_pos - r_pos_1), axes=2)
    # print(beta)
    beta = 1

    lv_vlm_img = imageData['lv_vlm_img']
    rv_vlm_img = imageData['rv_vlm_img']
    lv_p_img = imageData['lv_p_img']
    rv_p_img = imageData['rv_p_img']
    epi_p_img = imageData['epi_p_img']

    lv_vlm_cur = vals['lv_vlm_cur']
    rv_vlm_cur = vals['rv_vlm_cur']
    lv_p_cur = vals['lv_p_cur']
    rv_p_cur = vals['rv_p_cur']
    epi_p_cur = vals['epi_p_cur']

    print(vals['lv_vlm_dat'], vals['lv_vlm_ref'], lv_vlm_cur)
    print(vals['rv_vlm_dat'], vals['rv_vlm_ref'], rv_vlm_cur)
    print(lv_vlm_img)
    print(rv_vlm_img)
    lv_p_err = 0
    rv_p_err = 0
    epi_p_err = 0
    a = lv_p_img[0][0,:]
    for i in range(len(lv_p_img)):
        for j in range(len(lv_p_img[i])):
            lv_p_err = lv_p_err + np.linalg.norm(lv_p_img[i][j, :] - lv_p_cur[i][j, :])
            rv_p_err = rv_p_err + np.linalg.norm(rv_p_img[i][j, :] - rv_p_cur[i][j, :])
            epi_p_err = epi_p_err + np.linalg.norm(epi_p_img[i][j, :] - epi_p_cur[i][j, :])
    lv_vlm_err = 0
    rv_vlm_err = 0
    for i in range(len(lv_vlm_img)):
        lv_vlm_err = lv_vlm_err + np.abs(lv_vlm_img[i] - lv_vlm_cur[i])
        rv_vlm_err = rv_vlm_err + np.abs(rv_vlm_img[i] - rv_vlm_cur[i])
    lv_vlm_err = lv_vlm_err / lv_vlm_img[-1] * 10
    rv_vlm_err = rv_vlm_err / rv_vlm_img[-1] * 10

    fit_err = lv_vlm_err + rv_vlm_err + (epi_p_err + lv_p_err + rv_p_err) * 0.1

    print("LV p err: " + str(lv_p_err))
    print("RV p err: " + str(rv_p_err))
    print("EPI p err: " + str(epi_p_err))
    print("LV v err: " + str(lv_vlm_err))
    print("RV v err: " + str(rv_vlm_err))
    print("fit LVRV err: " + str(lv_vlm_err + rv_vlm_err))
    print("fit tot err: " + str(fit_err))




    ####################### Finish Calculation  #######################################################
    cal_vals = {
        'beta': beta,
        'r_pos': r_pos,
        'r_disp_max': r_disp_max,
        'fit_err': fit_err
    }

    return cal_vals
