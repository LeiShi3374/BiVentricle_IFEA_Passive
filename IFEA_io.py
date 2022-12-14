## This file contains the input and output codes for vtk files

import vtk
from vtkmodules.util import numpy_support as vtknp

def read_surf(file_name):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(file_name)
    reader.Update()
    output = reader.GetOutput()
    gnode_vtk_array = output.GetPointData().GetArray('GlobalNodeID')
    gnode = vtknp.vtk_to_numpy(gnode_vtk_array) - 1
    return output, gnode

def read_fiber(file_name):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(file_name)
    reader.Update()
    long_output = reader.GetOutput()
    fiblong_vtk_array = long_output.GetCellData().GetArray('FIB_DIR')
    # fiblong_vtk_array = long_output.GetPointData().GetArray('FIB_DIR')
    fiblong = vtknp.vtk_to_numpy(fiblong_vtk_array)
    return fiblong

def read_vol(file_name):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(file_name)
    reader.Update()
    output = reader.GetOutput()
    return output

def write_vol(output, file_name):
    mshwrite = vtk.vtkXMLUnstructuredGridWriter()
    mshwrite.SetInputData(output)
    mshwrite.SetFileName(file_name)
    mshwrite.Write()

def write_surf(output, file_name):
    mshwrite = vtk.vtkXMLPolyDataWriter()
    mshwrite.SetInputData(output)
    mshwrite.SetFileName(file_name)
    mshwrite.Write()

def write_fiber(fiberdata, file_name, output):
    cur_NoC = len(fiberdata)
    bff = vtk.vtkDoubleArray()
    bff.SetNumberOfComponents(3)
    bff.Allocate(cur_NoC)
    bff.SetNumberOfTuples(cur_NoC)
    bff.SetName("FIB_DIR")
    for i in range(cur_NoC):
        bff.SetTuple3(i, fiberdata[i, 0], fiberdata[i, 1], fiberdata[i, 2])
    output.GetCellData().AddArray(bff)
    # output.GetPointData().AddArray(bff)
    mshwrite = vtk.vtkXMLUnstructuredGridWriter()
    mshwrite.SetInputData(output)
    mshwrite.SetFileName(file_name)
    mshwrite.Write()