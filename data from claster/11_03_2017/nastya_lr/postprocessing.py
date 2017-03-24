import numpy as np
import pandas as pd
import vtk
from vtk import *
from vtk.util import numpy_support
import matplotlib.pylab as plt
import os
import imp
import h5py
from sklearn.neighbors import KDTree
from tqdm import tqdm

FILE_ECG = 'data_1.ecg'
FILE_CONFIG = 'config.py'
FILE_MESH = 'init_mesh_1.vtu'
FILE_H5 = 'results.h5'


def _check_int(s):
    if s[0] in ('-', '+'):
        return s[1:].isdigit()
    return s.isdigit()


class Maps:
    """
    Write activation, repolarization and APD90 maps in Chaste .vtu output
    """

    def __init__(self, file_in, file_out=None):
        """
        Constructor
        :param file_in: Input file (.vtu)
        :param file_out: Output file (.vtu)
        """
        self.file_in = file_in
        self.file_out = file_out
        self._apd_foot = None
        self._apd_overshoot = None
        self._apd_treshold = None
        self._set_apd_treshold(-86.03955628, 29.77557964)

        self.reader = None
        self.data_heart = None
        self.points = None
        self.cells = None
        self.arr_number = None
        self.activation_map = None
        self.repolarization_map = None
        self.time = None
        self.early_point = None
        self.late_point = None
        self.fl = None

    def _set_apd_treshold(self, apd_foot, apd_treshold):
        """
        Set foot and overshoot values for find APD90 threshold.
        :param apd_foot: foot level APD
        :param apd_treshold: Overshoot maximum
        :return:
        """
        self._apd_foot = apd_foot
        self._apd_overshoot = apd_treshold
        self._apd_treshold = self._apd_foot + np.abs(self._apd_overshoot - self._apd_foot) * 0.1

    def update(self):
        """
        Run treatment
        :return:
        """
        assert (os.path.exists(self.file_in))
        assert (self._apd_treshold is not None)
        assert (self._apd_overshoot is not None)
        assert (self._apd_treshold is not None)

        self.reader = vtk.vtkXMLUnstructuredGridReader()
        self.reader.SetFileName(self.file_in)
        self.reader.Update()

        self.data_heart = self.reader.GetOutput()
        self.points = numpy_support.vtk_to_numpy(self.data_heart.GetPoints().GetData())
        self.cells = numpy_support.vtk_to_numpy(self.data_heart.GetCells().GetData())

        self.fl = h5py.File(FILE_H5, 'r')
        self.fl.items()
        self.data = self.fl['Data']


        self.arr_number = self.data_heart.GetPointData().GetNumberOfArrays()
        self.activation_map = np.zeros(len(self.points))
        self.repolarization_map = np.zeros(len(self.points))
        self.activation_map[:] = np.nan
        self.repolarization_map[:] = np.nan
        self.time = 0
        print self.arr_number, len(self.points), len(self.data), len(self.data[0])
        #data[:, 100, 0]
        print self.data_heart.GetPointData()
        #for i in xrange(self.arr_number):
        #for i in xrange(len(self.data)):
        #    if 'V_' in self.data_heart.GetPointData().GetArrayName(i):
        #        #data = numpy_support.vtk_to_numpy(self.data_heart.GetPointData().GetArray(i))
        #        data = self.data[:, i, 0]
        #        for j in xrange(len(data)):
        #            if data[j] > self._apd_treshold and np.isnan(self.activation_map[j]) and np.isnan(
        #                    self.repolarization_map[j]):
        #                self.activation_map[j] = self.time
        #            if data[j] < self._apd_treshold and not np.isnan(self.activation_map[j]) and np.isnan(
        #                    self.repolarization_map[j]):
        #                self.repolarization_map[j] = self.time
        #        self.time += 1

        for j in tqdm(xrange(len(self.data))):
            for i in xrange(len(self.data[0])):
                for_point = self.data[:, i, 0]
                if for_point[self.time] > self._apd_treshold and np.isnan(self.activation_map[i]) and np.isnan(
                        self.repolarization_map[i]):
                    self.activation_map[i] = self.time
                if for_point[self.time] < self._apd_treshold and not np.isnan(self.activation_map[i]) and np.isnan(
                        self.repolarization_map[i]):
                    self.repolarization_map[i] = self.time
            self.time += 1

        self.early_point = np.argmin(self.activation_map)
        self.late_point = np.argmax(self.activation_map)

        if self.file_out is not None and os.path.exists(os.path.dirname(self.file_out)):
            array_for_write = numpy_support.numpy_to_vtk(self.activation_map, deep=True, array_type=vtk.VTK_FLOAT)
            array_for_write.SetName("activation_map")
            self.data_heart.GetPointData().AddArray(array_for_write)

            array_for_write2 = numpy_support.numpy_to_vtk(self.repolarization_map, deep=True, array_type=vtk.VTK_FLOAT)
            array_for_write2.SetName("repolarization_map")
            self.data_heart.GetPointData().AddArray(array_for_write2)

            array_for_write3 = numpy_support.numpy_to_vtk(self.repolarization_map - self.activation_map, deep=True,
                                                          array_type=vtk.VTK_FLOAT)
            array_for_write3.SetName("APD90_map")
            self.data_heart.GetPointData().AddArray(array_for_write3)

            writer = vtk.vtkXMLUnstructuredGridWriter()
            writer.SetFileName(self.file_out)
            writer.SetInputData(self.data_heart)

            writer.Write()

    def report(self, plots=False):
        """
        Write summary information about simulation data
        :param plots: In case of True, showed plots
        :return:
        """
        print "Activation map (array, min, max)", self.activation_map, min(self.activation_map), max(
            self.activation_map)
        print "Repolarization map (array, min, max)", self.repolarization_map, min(self.repolarization_map), max(
            self.repolarization_map)
        print "APD90 map (array, min, max, rep. disp.)", min(self.repolarization_map), max(
            self.repolarization_map), max(self.repolarization_map) - min(self.repolarization_map)
        print "Early and late points", self.early_point, self.late_point
        if plots:
            plt.hist(self.activation_map, bins=30, normed=True)
            plt.title("Activation")
            plt.show
            plt.close()
            plt.hist(self.repolarization_map, bins=30, normed=True)
            plt.title("Repolarization")
            plt.show()
            plt.close()
            plt.hist(self.repolarization_map - self.activation_map, bins=30, normed=True)
            plt.title("APD90")
            plt.show()
            plt.close()


class ModelECG:
    """
    Write files with ECG
    """

    def __init__(self, folder_data, folder_in):
        self.folder_data = folder_data
        self.config_file = 'config.py'
        self.folder_in = folder_in

    def update(self):
        assert (os.path.exists(os.path.join(self.folder_data, self.config_file)))
        assert (os.path.join(self.folder_in, FILE_MESH))

        # Load config file for leads placement
        config = imp.load_source('*', os.path.join(self.folder_data, self.config_file))

        # Open mesh
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(os.path.join(self.folder_in, FILE_MESH))
        reader.Update()

        data_heart = reader.GetOutput()
        points = numpy_support.vtk_to_numpy(data_heart.GetPoints().GetData())
        cells = numpy_support.vtk_to_numpy(data_heart.GetCells().GetData())

        assert (len(points) > 0)
        assert (len(cells) > 0)

        torso_tree = KDTree(points, leaf_size=2)
        fl = h5py.File(os.path.join(self.folder_in, FILE_H5), 'r')

        temp = torso_tree.query_radius(np.array(config.lead_left).reshape(1, -1), r=config.radius)[0]
        self.potential_left = -fl['Data'][:, temp[0], 1]
        temp = torso_tree.query_radius(np.array(config.lead_right).reshape(1, -1), r=config.radius)[0]
        self.potential_right = -fl['Data'][:, temp[0], 1]
        temp = torso_tree.query_radius(np.array(config.lead_foot).reshape(1, -1), r=config.radius)[0]
        self.potential_foot = -fl['Data'][:, temp[0], 1]

        self.ecg_I = (self.potential_right - self.potential_left) / 2
        self.ecg_II = (self.potential_right - self.potential_foot) / 2.
        self.ecg_III = (self.potential_left - self.potential_foot) / 2.

    def get_frontal(self):
        assert (self.ecg_I is not None)
        return pd.DataFrame({"I": self.ecg_I,
                             "II": self.ecg_II,
                             "III": self.ecg_III})


class AmikardECG:
    """
    Load and work with Amikard ecg format
    """

    def __init__(self, file):
        self.file = file
        self.ecg = None

    def update(self):
        dictionary = {}
        key = None

        with open(self.file, 'r') as fl:
            for line in fl:
                if '*ECG' in line:
                    key = line.split(' ')[1]
                    dictionary[key] = []
                elif _check_int(line.strip()):
                    dictionary[key].append(int(line.strip()) / 1000.)
                else:
                    pass

        self.ecg = pd.DataFrame(dictionary)

    def get_frontal(self):
        assert (self.ecg is not None)
        return self.ecg


class RealDataECGComparator:
    """
    Compare amikard ecg data with simulation
    """

    def __init__(self, real_ecg, model_ecg, shift=50, folder_out=None):
        self.real_ecg = real_ecg
        self.model_ecg = model_ecg
        self.shift = shift
        self.folder_out = folder_out
        self.ax_I = None
        self.ax_II = None
        self.ax_III = None
        self.filename_I = 'lead_I.png'
        self.filename_II = 'lead_II.png'
        self.filename_III = 'lead_III.png'

    def update(self):
        assert (type(self.real_ecg) == pd.DataFrame)
        assert (type(self.model_ecg) == pd.DataFrame)

        self.f_I, self.ax_I = plt.subplots()
        self.ax_I.plot(self.real_ecg['I'], color='black')
        self.ax_I.plot(np.append(np.zeros(self.shift), self.model_ecg['I']), color='red')
        self.ax_I.set_ylim((-2, 2))

        self.f_II, self.ax_II = plt.subplots()
        self.ax_II.plot(self.real_ecg['II'], color='black')
        self.ax_II.plot(np.append(np.zeros(self.shift), self.model_ecg['II']), color='green')
        self.ax_II.set_ylim((-2, 2))

        self.f_III, self.ax_III = plt.subplots()
        self.ax_III.plot(self.real_ecg['III'], color='black')
        self.ax_III.plot(np.append(np.zeros(self.shift), self.model_ecg['III']), color='blue')
        self.ax_III.set_ylim((-2, 2))

        if self.folder_out is not None:
            self.f_I.savefig(os.path.join(self.folder_out, self.filename_I))
            self.f_II.savefig(os.path.join(self.folder_out, self.filename_II))
            self.f_III.savefig(os.path.join(self.folder_out, self.filename_III))


class ECG:
    def __init__(self, folder_config, folder_in, folder_out):
        self.folder_config = folder_config
        self.folder_in = folder_in
        self.folder_out = folder_out
        self.shift = 150

    def update(self):
        assert (os.path.exists(self.folder_in))
        assert (os.path.exists(self.folder_out))
        assert (os.path.exists(self.folder_config))

        model_ECG = ModelECG(self.folder_config, self.folder_in)
        model_ECG.update()
        model_frontal = model_ECG.get_frontal()

        real_ECG = AmikardECG(os.path.join(self.folder_config, FILE_ECG))
        real_ECG.update()
        real_frontal = real_ECG.get_frontal()

        cmp = RealDataECGComparator(real_frontal, model_frontal, shift=self.shift, folder_out=self.folder_out)
        cmp.update()
