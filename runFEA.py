'''This file is to update the svFSI.inp file and run the simulation
 to generate the simulation result
'''

import numpy as np
import os
import glob

def runFEA(mat_para_opt, f_ind, mat_para_fix, finalflag):
    ######## Generate New svFSI file ##############################
    fname = os.getcwd() + "/"
    m_a = mat_para_opt[0]
    m_b = mat_para_opt[1]
    m_af = mat_para_opt[2]
    m_bf = mat_para_opt[3]
    m_as = mat_para_opt[4]
    m_bs = mat_para_opt[5]
    m_afs = mat_para_opt[6]
    m_bfs = mat_para_opt[7]

    m_eta = mat_para_fix[0]
    m_E = mat_para_fix[1]
    # lv_load = mat_para_fix[2]
    # rv_load = mat_para_fix[3]

    k = '%02d' % (f_ind)
    in_name = "US_svFSI_" + "ref.inp"
    out_name = "US_svFSI_" + k + ".inp"

    with open(in_name, 'r') as svFile:
        svRead = svFile.readlines()
    for item in range(len(svRead)):
        # if "Constitutive model" in svRead[item]:
        if "Elasticity modulus: " in svRead[item]:
            svRead[item] = "   Elasticity modulus: " + str(m_E) + "\n"
        if "Viscosity: " in svRead[item]:
            svRead[item] = "   Viscosity: " + str(m_eta) + "\n"
        if " a: " in svRead[item]:
            svRead[item] = "      a: " + str(m_a) + "\n"
        if " b: " in svRead[item]:
            svRead[item] = "      b: " + str(m_b) + "\n"
        if " a4f: " in svRead[item]:
            svRead[item] = "      a4f: " + str(m_af) + "\n"
        if " b4f: " in svRead[item]:
            svRead[item] = "      b4f: " + str(m_bf) + "\n"
        if " a4s: " in svRead[item]:
            svRead[item] = "      a4s: " + str(m_as) + "\n"
        if " b4s: " in svRead[item]:
            svRead[item] = "      b4s: " + str(m_bs) + "\n"
        if " afs: " in svRead[item]:
            svRead[item] = "      afs: " + str(m_afs) + "\n"
        if " bfs: " in svRead[item]:
            svRead[item] = "      bfs: " + str(m_bfs) + "\n"

    for item in range(len(svRead)):
        if "Face file path: " in svRead[item]:
            svReadsplt = svRead[item].split("/")
            svRead[item] = svRead[item].replace(svReadsplt[0], "      Face file path: " + "05_mesh_" + k)
        if "Mesh file path: " in svRead[item]:
            svReadsplt = svRead[item].split("/")
            svRead[item] = svRead[item].replace(svReadsplt[0], "      Mesh file path: " + "05_mesh_" + k)
        if "Fiber direction file path: " in svRead[item]:
            svReadsplt = svRead[item].split("/")
            svRead[item] = svRead[item].replace(svReadsplt[0], "      Fiber direction file path: " + "05_mesh_" + k)
        if "Save results in folder: " in svRead[item]:
            svRead[item] = "Save results in folder: " + "05_mesh_" + k + "result\n"

        if finalflag == True:
            if "Continue previous simulation: " in svRead[item]:
                svRead[item] = "Continue previous simulation: 1\n"
            if "Number of time steps" in svRead[item]:
                svRead[item] = "Number of time steps: 20\n"

    with open(out_name, 'w') as svFileNew:
        svFileNew.writelines(svRead)

    # with open("load.dat", "r") as load_file:
    #     load_read = load_file.readlines()
    #     load_read[2] = "0.43" + "  " + str(lv_load)

    # with open("load_r.dat", "r") as load_rfile:
    #     load_rread = load_rfile.readlines()
    #     load_rread[2] = "0.43" + "  " + str(rv_load)

    # with open("load.dat", 'w') as load_fileNew:
    #     load_fileNew.writelines(load_read)

    # with open("load_r.dat", 'w') as load_rfileNew:
    #     load_rfileNew.writelines(load_rread)

    if finalflag == False:
        rmfd = glob.glob("05_mesh_" + k + "result/*")
        for f in rmfd:
            os.remove(f)

    os.system("mpiexec -np 8 '/home/lei/svFSI_build_20220823/svFSI-build/bin/svFSI' " + out_name)

    ######## Finish Generate New svFSI file ##############################

    return 0
