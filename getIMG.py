## This is used to calculate the fluid volume encompassed by a surface
## Currently a divergence theory is used

import numpy as np
import vtk
from vtkmodules.util import numpy_support as vtknp
from FluidVolume import FluidVolume


def getIMG(imagePath, file_index, vals, p_info):
    # file_index = ["01_RR_70", "02_RR_80", "03_RR_90", "04_RR_99", "05_RR_10", "06_RR_20", "07_RR_30", "08_RR_40",
    #               "09_RR_50",
    #               "10_RR_60"]

    lv_p_gindex = p_info[0]
    rv_p_gindex = p_info[1]
    epi_p_gindex = p_info[2]

    lv_vlm_img = np.zeros(len(file_index))
    rv_vlm_img = np.zeros(len(file_index))
    lv_p_img = [None] * len(file_index)
    rv_p_img = [None] * len(file_index)
    epi_p_img = [None] * len(file_index)

    for f_index in range(len(file_index)):
        # LV Volume
        file_name = imagePath + "endo_lv/" + file_index[f_index] + "_registered.vtk"
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(file_name)
        reader.Update()
        lv_i_output = reader.GetOutput()
        point0img = vtknp.vtk_to_numpy(lv_i_output.GetPoints().GetData())
        lv_gnode = vals['lv_gnode']
        lv_edge_gnode = vals['lv_edge_gnode']
        lv_edge_lnode = np.zeros(len(lv_edge_gnode)).astype(int)
        for i in range(len(lv_edge_gnode)):
            lv_edge_lnode[i] = np.where(lv_gnode == lv_edge_gnode[i])[0]
        lv_i_NoC = lv_i_output.GetNumberOfCells()
        conn_i_lv = []
        for i in range(lv_i_NoC):
            cell = lv_i_output.GetCell(i)
            npts = cell.GetNumberOfPoints()
            connt = [cell.GetPointId(j) for j in range(npts)]
            conn_i_lv.append(connt)
        conn_i_lv = np.asarray(conn_i_lv)
        lv_img_lnode = np.linspace(0, len(point0img) - 1, num=len(point0img)).astype(int)

        lv_vlm_img[f_index] = FluidVolume(point0img, lv_edge_lnode, lv_img_lnode, conn_i_lv)  # lv_vlm_img

        lv_p_img[f_index] = np.zeros((len(lv_p_gindex), 3))
        lv_p_lindex = np.zeros(len(lv_p_gindex)).astype(int)
        for i in range(len(lv_p_gindex)):
            lv_p_lindex[i] = np.where(lv_gnode == lv_p_gindex[i])[0]
            lv_p_img[f_index][i, :] = point0img[lv_p_lindex[i], :]

        # RV Volume
        file_name = imagePath + "endo_rv/" + file_index[f_index] + "_registered.vtk"
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(file_name)
        reader.Update()
        rv_i_output = reader.GetOutput()
        point0img = vtknp.vtk_to_numpy(rv_i_output.GetPoints().GetData())
        rv_gnode = vals['rv_gnode']
        rv_edge_gnode = vals['rv_edge_gnode']
        rv_edge_lnode = np.zeros(len(rv_edge_gnode)).astype(int)
        for i in range(len(rv_edge_gnode)):
            rv_edge_lnode[i] = np.where(rv_gnode == rv_edge_gnode[i])[0]
        rv_i_NoC = rv_i_output.GetNumberOfCells()
        conn_i_rv = []
        for i in range(rv_i_NoC):
            cell = rv_i_output.GetCell(i)
            npts = cell.GetNumberOfPoints()
            connt = [cell.GetPointId(j) for j in range(npts)]
            conn_i_rv.append(connt)
        conn_i_rv = np.asarray(conn_i_rv)
        rv_img_lnode = np.linspace(0, len(point0img) - 1, num=len(point0img)).astype(int)

        rv_vlm_img[f_index] = FluidVolume(point0img, rv_edge_lnode, rv_img_lnode, conn_i_rv)  # rv_vlm_img

        rv_p_img[f_index] = np.zeros((len(rv_p_gindex), 3))
        rv_p_lindex = np.zeros(len(rv_p_gindex)).astype(int)
        for i in range(len(rv_p_gindex)):
            rv_p_lindex[i] = np.where(rv_gnode == rv_p_gindex[i])[0]
            rv_p_img[f_index][i, :] = point0img[rv_p_lindex[i], :]

        # epi
        file_name = imagePath + "epi/" + file_index[f_index] + "_registered.vtk"
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(file_name)
        reader.Update()
        epi_i_output = reader.GetOutput()
        point0img = vtknp.vtk_to_numpy(epi_i_output.GetPoints().GetData())
        epi_gnode = vals['epi_gnode']

        epi_p_img[f_index] = np.zeros((len(epi_p_gindex), 3))
        epi_p_lindex = np.zeros(len(epi_p_gindex)).astype(int)
        for i in range(len(epi_p_gindex)):
            epi_p_lindex[i] = np.where(epi_gnode == epi_p_gindex[i])[0]
            epi_p_img[f_index][i, :] = point0img[epi_p_lindex[i], :]

    imageData = {
        'lv_vlm_img': lv_vlm_img,
        'rv_vlm_img': rv_vlm_img,
        'rv_p_img': rv_p_img,
        'lv_p_img': lv_p_img,
        'epi_p_img': epi_p_img,
    }

    return imageData
