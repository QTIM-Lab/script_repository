import glob
import fnmatch
import os
import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib
from shutil import copy, move
from qtim_tools.qtim_utilities import nifti_util
# import nifti_util
import csv

def Repeatability_Assay():

	CED_Patients = ['CED_' + str(item).zfill(6) for item in np.arange(1,43)]
	NHX_Patients = ['NHX_' + str(item).zfill(6) for item in np.arange(1,9)]
	ALL_Patients = CED_Patients + NHX_Patients

	blur_params = ['blur_' + str(item).zfill(6) for item in [0,0.2,0.8,1.2]]
	pca_params = ['pca_' + str(item).zfill(6) for item in [0,1,2,3,4]]
	noise_threshold_params = ['threshold_' + str(item).zfill(6) for item in [-1, 0.01]]
	integration_params = ['recursive_', 'conv_']

	# on/off params
	aif_params = ['autoAIF', '']
	t1map_params = ['t1map', '']

	visits = ['VISIT_01', 'VISIT_02']

	suffix = 'Function_Test_ktrans.nii.gz'

	for patient in ALL_Patients:

		for aif in aif_params:
			for t1map in t1map_params:
				for integration in integration_params:
					for blur in blur_params:
						for pca in pca_params:
							for threshold in noise_threshold_params:

								filename = '_'.join(patient, visit, aif, t1map, integration, blur, pca, threshold, suffix)
								

								try:

								except:
									'No such file at ' + filename

	return

if __name__ == '__main__':
	Repeatability_Assay()