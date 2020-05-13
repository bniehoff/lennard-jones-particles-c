# Microcanonical Ensemble Simulation
# batch-simulate.py:  Batch simulation script
# 
# (c) 2020 Ben Niehoff
# 
# Description:
# This script runs the main program md-simulate over a spread of
# parameters in order to generate data which can then be examined in
# the Mathematica notebook Visualization.nb

import os
import sys
import shutil
import getopt
import ast
import numpy as np
import pandas as pd

# Get options
opts, args = getopt.getopt(sys.argv[1:], 'rc:d:t:o:',
	['remove-data', 'cellcount=', 'density-linspace=', 'temperature-linspace=',
	'output-directory='])

remove_data = False
cellcount = 5
density_linspace = [0.8, 0.8, 1]
# temperature_linspace = [0.9, 0.9, 1]
temperature_linspace = [0.1, 0.9, 5]
# temperature_linspace = [0.1, 0.9, 17]
# temperature_linspace = [0.1, 0.9, 33]
# temperature_linspace = [0.125, 0.875, 16]
output_directory = 'data'

#print(opts)

for opt, val in opts:
	if opt in ['-r', '--remove-data']:
		remove_data = True
	elif opt in ['-c', '--cellcount']:
		cellcount = int(val)
	elif opt in ['-d', '--density-linspace']:
		density_linspace = ast.literal_eval(val)
	elif opt in ['-t', '--temperature-linspace']:
		temperature_linspace = ast.literal_eval(val)
	elif opt in ['-o', '--output-directory']:
		output_directory = val

# First make sure ./data directory exists
if not os.path.exists(output_directory):
	print('Batch simulate:', output_directory, 'does not exist')
	print('Batch simulate: Creating', output_directory)
	if remove_data:
		print('Batch simulate: No need to remove files from', output_directory)
		remove_data = False
	os.mkdir(output_directory)

# Check whether to clean the data directory or not
if remove_data:
	# Clean output directory of all files
	print('Batch simulate: Removing files from', output_directory)
	for root, dirs, files in os.walk(output_directory):
		for f in files:
			os.remove(os.path.join(root, f))
		for d in dirs:
			shutil.rmtree(os.path.join(root, d))

# Prepare initial file
if not os.path.exists(os.path.join(output_directory, 'thermo_measurements.csv')):
	thermo_meas_file = open(os.path.join(output_directory, 'thermo_measurements.csv'), 'w')
	thermo_meas_file.write('Density,Temp,Energy,HeatCapCv,Pressure\n')
	thermo_meas_file.close()

# Run the code
print('Batch simulate: Ready to run simulation for densities '
	'in {} and temperatures in {}'.format(density_linspace, temperature_linspace))
print('Batch simulate: Beginning simulation now...')

#exit()

densities = np.linspace(*density_linspace)
temperatures = np.linspace(*temperature_linspace)

#print(densities)
#print(temperatures

#temperatures = np.linspace(0.1,0.9,num=17)
#temperatures = [0.1, 0.3, 0.5, 0.7, 0.9]

# Run the simulations
for rho in densities:
	rho_dirname = 'rho_{:.3f}'.format(rho)
	if not os.path.exists(os.path.join(output_directory, rho_dirname)):
		os.mkdir(os.path.join(output_directory, rho_dirname))
	
	for T in temperatures:
		T_dirname = 'T_{:.3f}'.format(T)
		
		if not os.path.exists(os.path.join(output_directory, rho_dirname, T_dirname)):
			os.mkdir(os.path.join(output_directory, rho_dirname, T_dirname))
		
		if not os.path.exists(os.path.join(output_directory, rho_dirname, T_dirname, 'time_series.csv')):
			time_series_file = open(os.path.join(output_directory, rho_dirname, T_dirname, 'time_series.csv'), 'w')
			time_series_file.write('TimeStep,Temp,PotEnergy,TotEnergy,MeanSqDisp\n')
			time_series_file.close()
		
		if not os.path.exists(os.path.join(output_directory, rho_dirname, T_dirname, 'final_state.csv')):
			final_state_file = open(os.path.join(output_directory, rho_dirname, T_dirname, 'final_state.csv'), 'w')
			final_state_file.write('PosX,PosY,PosZ,VelX,VelY,VelZ,Speed\n')
			final_state_file.close()
		
		if not os.path.exists(os.path.join(output_directory, rho_dirname, T_dirname, 'summary_info.csv')):
			summary_info_file = open(os.path.join(output_directory, rho_dirname, T_dirname, 'summary_info.csv'), 'w')
			summary_info_file.write('CellCount,SideLength,AtomCount,BlockCount,BlocksPerSide,BlockLength\n')
			summary_info_file.close()
		
		cmd = ' '.join(['./md-simulate',
			'--cellcount', str(cellcount),
			'--density', str(rho),
			'--temperature', str(T),
			'--output-directory', output_directory,
			'--prefix', '"' + rho_dirname + '/' + T_dirname + '"'])
		
		print('Batch simulate:  Executing with density {:.3f} and temperature {:.3f}'.format(rho, T))
		#print(cmd)
		os.system(cmd)

# os.system('./md-simulate --cellcount 5 --density 0.8 --temperature 0.9 --output-directory data --prefix T0.9')
