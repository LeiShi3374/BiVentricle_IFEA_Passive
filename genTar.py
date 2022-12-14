## This file gets the vtk structures and values from the calculation
## to generate the new FEA files

import numpy as np
import os
import vtk
from vtkmodules.util import numpy_support as vtknp
import IFEA_io as Iio

def genTar(fname, f_ind, vals, cal_vals):

    point0d = vals['point0d']
    point0 = vals['point0']
    lv_output = vals['lv_output']
    lv_gnode = vals['lv_gnode']
    rv_output = vals['rv_output']
    rv_gnode = vals['rv_gnode']
    top_output = vals['top_output']
    top_gnode = vals['top_gnode']
    epi_output = vals['epi_output']
    epi_gnode = vals['epi_gnode']
    epi_apex_output = vals['epi_apex_output']
    epi_apex_gnode = vals['epi_apex_gnode']
    epi_mid_output = vals['epi_mid_output']
    epi_mid_gnode = vals['epi_mid_gnode']
    ref_output = vals['ref_output']
    cur_output = vals['cur_output']

    beta = cal_vals['beta']
    r_pos = cal_vals['r_pos']

    ################## Read and write vtu files ######################################################
    ## write vtu files
    # disp_abs = point1 - point0d

    k = '%02d' % (f_ind + 1)
    fname_ref = "05_mesh_ref/"
    fname_tar = "05_mesh_" + k + "/"
    file_target = fname + fname_tar + "mesh-complete.mesh.vtu"
    file_target_lv = fname + fname_tar + "mesh-surfaces/endo_lv.vtp"
    file_target_rv = fname + fname_tar + "mesh-surfaces/endo_rv.vtp"
    file_target_top = fname + fname_tar + "mesh-surfaces/top.vtp"
    file_target_long = fname + fname_tar + "fibersLong.vtu"
    file_target_sheet = fname + fname_tar + "fibersSheet.vtu"
    file_target_epi = fname + fname_tar + "mesh-surfaces/epi.vtp"
    file_target_epi_apex = fname + fname_tar + "mesh-surfaces/epi_apex.vtp"
    file_target_epi_mid = fname + fname_tar + "mesh-surfaces/epi_mid.vtp"

    file_long = fname + fname_ref + "fibersLong.vtu"
    file_sheet = fname + fname_ref + "fibersSheet.vtu"

    if os.path.exists(fname + fname_tar) == False:
        os.mkdir(fname + fname_tar)
    if os.path.exists(fname + fname_tar + "mesh-surfaces/") == False:
        os.mkdir(fname + fname_tar + "mesh-surfaces/")

    point2 = point0 - beta * r_pos
    point2_vtk = vtk.vtkPoints()

    F_vtk_array = cur_output.GetPointData().GetArray('Def_grad')
    F = vtknp.vtk_to_numpy(F_vtk_array)

    fiblong = Iio.read_fiber(file_long)
    fibsheet = Iio.read_fiber(file_sheet)

    cur_NoC = cur_output.GetNumberOfCells()
    cur_NoP = cur_output.GetNumberOfPoints()
    conn_cur = []
    for i in range(cur_NoC):
        cell = cur_output.GetCell(i)
        npts = cell.GetNumberOfPoints()
        connt = [cell.GetPointId(j) for j in range(npts)]
        conn_cur.append(connt)
    conn_cur = np.asarray(conn_cur)

    fiblong2 = np.zeros((cur_NoC, 3))
    fibsheet2 = np.zeros((cur_NoC, 3))
    for i in range(cur_NoC):
        PA0 = point0d[conn_cur[i, 0], :]
        PB0 = point0d[conn_cur[i, 1], :]
        PC0 = point0d[conn_cur[i, 2], :]
        PD0 = point0d[conn_cur[i, 3], :]
        PA2 = point2[conn_cur[i, 0], :]
        PB2 = point2[conn_cur[i, 1], :]
        PC2 = point2[conn_cur[i, 2], :]
        PD2 = point2[conn_cur[i, 3], :]
        FVF0 = fiblong[i, :]
        FVS0 = fibsheet[i, :]
        MP0 = np.array((PB0 - PA0, PC0 - PA0, PD0 - PA0)).T
        NF0 = np.linalg.solve(MP0, FVF0)
        NS0 = np.linalg.solve(MP0, FVS0)
        fiblong2[i, :] = (PB2 - PA2) * NF0[0] + (PC2 - PA2) * NF0[1] + (PD2 - PA2) * NF0[2]
        fiblong2[i, :] = fiblong2[i, :] / np.linalg.norm(fiblong2[i, :])
        fibsheet2[i, :] = (PB2 - PA2) * NS0[0] + (PC2 - PA2) * NS0[1] + (PD2 - PA2) * NS0[2]
        fibsheet2[i, :] = fibsheet2[i, :] / np.linalg.norm(fibsheet2[i, :])

    # fiblong2 = np.zeros((cur_NoP, 3))
    # fibsheet2 = np.zeros((cur_NoP, 3))
    # for i in range(cur_NoP):
    #     fiblong2[i, 0] = F[i, 0] * fiblong[i, 0] + F[i, 1] * fiblong[i, 1] + F[i, 2] * fiblong[i, 2]
    #     fiblong2[i, 1] = F[i, 3] * fiblong[i, 0] + F[i, 4] * fiblong[i, 1] + F[i, 5] * fiblong[i, 2]
    #     fiblong2[i, 2] = F[i, 6] * fiblong[i, 0] + F[i, 7] * fiblong[i, 1] + F[i, 8] * fiblong[i, 2]
    #     fiblong2[i, :] = fiblong2[i, :] / np.linalg.norm(fiblong2[i, :])
    #     fibsheet2[i, 0] = F[i, 0] * fibsheet[i, 0] + F[i, 1] * fibsheet[i, 1] + F[i, 2] * fibsheet[i, 2]
    #     fibsheet2[i, 1] = F[i, 3] * fibsheet[i, 0] + F[i, 4] * fibsheet[i, 1] + F[i, 5] * fibsheet[i, 2]
    #     fibsheet2[i, 2] = F[i, 6] * fibsheet[i, 0] + F[i, 7] * fibsheet[i, 1] + F[i, 8] * fibsheet[i, 2]
    #     fibsheet2[i, :] = fibsheet2[i, :] / np.linalg.norm(fibsheet2[i, :])

    for i in range(len(point2)):
        point2_vtk.InsertNextPoint(point2[i, 0], point2[i, 1], point2[i, 2])
    ref_output.SetPoints(point2_vtk)

    point2_lv_vtk = vtk.vtkPoints()
    for i in range(len(lv_gnode)):
        point2_lv_vtk.InsertNextPoint(point2[lv_gnode[i], 0], point2[lv_gnode[i], 1], point2[lv_gnode[i], 2])
    lv_output.SetPoints(point2_lv_vtk)

    point2_rv_vtk = vtk.vtkPoints()
    for i in range(len(rv_gnode)):
        point2_rv_vtk.InsertNextPoint(point2[rv_gnode[i], 0], point2[rv_gnode[i], 1], point2[rv_gnode[i], 2])
    rv_output.SetPoints(point2_rv_vtk)

    point2_top_vtk = vtk.vtkPoints()
    for i in range(len(top_gnode)):
        point2_top_vtk.InsertNextPoint(point2[top_gnode[i], 0], point2[top_gnode[i], 1], point2[top_gnode[i], 2])
    top_output.SetPoints(point2_top_vtk)

    point2_epi_vtk = vtk.vtkPoints()
    for i in range(len(epi_gnode)):
        point2_epi_vtk.InsertNextPoint(point2[epi_gnode[i], 0], point2[epi_gnode[i], 1], point2[epi_gnode[i], 2])
    epi_output.SetPoints(point2_epi_vtk)

    point2_epi_apex_vtk = vtk.vtkPoints()
    for i in range(len(epi_apex_gnode)):
        point2_epi_apex_vtk.InsertNextPoint(point2[epi_apex_gnode[i], 0], point2[epi_apex_gnode[i], 1],
                                            point2[epi_apex_gnode[i], 2])
    epi_apex_output.SetPoints(point2_epi_apex_vtk)

    point2_epi_mid_vtk = vtk.vtkPoints()
    for i in range(len(epi_mid_gnode)):
        point2_epi_mid_vtk.InsertNextPoint(point2[epi_mid_gnode[i], 0], point2[epi_mid_gnode[i], 1],
                                           point2[epi_mid_gnode[i], 2])
    epi_mid_output.SetPoints(point2_epi_mid_vtk)

    ## Write updated data into files
    Iio.write_vol(ref_output, file_target)
    Iio.write_surf(lv_output, file_target_lv)
    Iio.write_surf(rv_output, file_target_rv)
    Iio.write_surf(top_output, file_target_top)
    Iio.write_surf(epi_output, file_target_epi)
    Iio.write_surf(epi_apex_output, file_target_epi_apex)
    Iio.write_surf(epi_mid_output, file_target_epi_mid)
    Iio.write_fiber(fiblong2, file_target_long, ref_output)
    Iio.write_fiber(fibsheet2, file_target_sheet, ref_output)

    ################## Finish Read and write vtu files ###############################################

    return 0