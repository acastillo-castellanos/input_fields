#!/usr/bin/env python

import numpy as np
import h5py

file_path = 'snapshot_t2.00000.h5'

def output_matrix(f, n, x, y, file):
    with open(file, 'wb') as fileID:
        # Write 'n' as double precision (float64)
        np.array([n], dtype=np.float64).tofile(fileID)
        
        # Write 'y' values as double precision
        for j in range(n):
            np.array([y[j]], dtype=np.float64).tofile(fileID)
        
        # Write 'x' values and the corresponding 'f' matrix values as double precision
        for i in range(n):
            np.array([x[i]], dtype=np.float64).tofile(fileID)
            for j in range(n):
                np.array([f[j, i]], dtype=np.float64).tofile(fileID)

# Open the HDF5 file
with h5py.File(file_path, 'r') as hdf:
    # Access groups and datasets
    cells_group = hdf['/Cells']
    geometry_group = hdf['/Geometry']
    topology_dataset = hdf['/Topology']

    # Access datasets within groups
    cells_f_dataset = cells_group['f']
    cells_p_dataset = cells_group['p']
    cells_ux_dataset = cells_group['u.x']
    geometry_points_dataset = geometry_group['Points']

    # Read the data from datasets
    cells_f_data = cells_f_dataset[:]
    cells_p_data = cells_p_dataset[:]
    cells_ux_data = cells_ux_dataset[:]
    geometry_points_data = geometry_points_dataset[:]
    topology_data = topology_dataset[:]

    # Print the shape of each dataset
    print("Cells/f dataset shape:", cells_f_data.shape)
    print("Cells/p dataset shape:", cells_p_data.shape)
    print("Cells/u.x dataset shape:", cells_ux_data.shape)
    print("Geometry/Points dataset shape:", geometry_points_data.shape)
    print("Topology dataset shape:", topology_data.shape)

from scipy.interpolate import griddata
cell_centers = np.mean(geometry_points_data[topology_data,:], axis=1)
cells_x, cells_y, cells_z = np.split(cell_centers,3,axis=1)

# Define a regular grid
L0 = 0.24
n = 1024
xc = ((np.arange(n)+0.5)/n - 0.5) * L0
yc = ((np.arange(n)+0.5)/n - 0.5) * L0
grid_x, grid_y = np.meshgrid(xc,yc)

# Perform the interpolation
grid_f = np.squeeze(griddata(cell_centers[:,0:2], cells_f_data, (grid_x, grid_y), method='linear', fill_value=0))
grid_p = np.squeeze(griddata(cell_centers[:,0:2], cells_p_data, (grid_x, grid_y), method='linear', fill_value=0))
grid_u = np.squeeze(griddata(cell_centers[:,0:2], cells_ux_data[:,0], (grid_x, grid_y), method='linear', fill_value=0))
grid_v = np.squeeze(griddata(cell_centers[:,0:2], cells_ux_data[:,1], (grid_x, grid_y), method='linear', fill_value=0))

output_matrix(grid_f, n, xc, yc, file_path[:-3]+'_f.bin')
output_matrix(grid_p, n, xc, yc, file_path[:-3]+'_p.bin')
output_matrix(grid_u, n, xc, yc, file_path[:-3]+'_u.bin')
output_matrix(grid_v, n, xc, yc, file_path[:-3]+'_v.bin')
