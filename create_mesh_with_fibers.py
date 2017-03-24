# coding=utf-8
__author__ = 'Anastasia Bazhutina'

from vtk import *
from vtk.util.numpy_support import vtk_to_numpy
from sklearn.neighbors import KDTree
from const import *
import numpy as np
import logging
import os




#функция создает и записывает .axi файл, в котором для каждой тетраэдральной ячейки указан вектор, лежащий в центре ячейки
def create_mesh_with_fibers():
    logger = logging.getLogger('simple_log')
    #чтение данных vtk с сеткой, создание октанта
    file_name = MESH_FILE_NAME
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(file_name)
    reader.Update()  # Needed because of GetScalarRange
    ug = reader.GetOutput()
    output = reader.GetOutput()
    print output
    boneLocator = vtk.vtkCellLocator()
    boneLocator.SetDataSet(output)
    boneLocator.BuildLocator()
    cell = vtk_to_numpy(ug.GetCells().GetData())
    print cell
    print len(cell)/4.

    #чтение данных ДТМРТ
    reader2 = vtkUnstructuredGridReader()
    way = FILE_DATA_VTK
    reader2.SetFileName(way)
    reader2.Update()
    ug2 = reader2.GetOutput()
    points2 = ug2.GetPoints()
    coords_vectors_data = vtk_to_numpy(ug2.GetPointData().GetArray('DTMRIFiberOrientation'))
    coord_points_data = vtk_to_numpy(points2.GetData())
    #построение дерева
    tree = KDTree(coord_points_data, leaf_size=5, metric='euclidean')



if __name__ == "__main__":
    create_mesh_with_fibers()