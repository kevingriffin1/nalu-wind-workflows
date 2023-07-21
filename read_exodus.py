import netCDF4
import numpy as np
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt

import os
os.chdir('/projects/hfm/kgriffin/nalu-wind-workflows/')

dir_names = ['dycoarse']; exo_filenames = ['bump.e']; nx = 212

ps = ['-r',':b']
pc = ['r','b'] 

fig, ax = plt.subplots(figsize=(8, 4))
for i_dir in range(len(dir_names)): 
    nc = netCDF4.Dataset('./'+dir_names[i_dir]+'/results/'+exo_filenames[i_dir])
    x_1d = nc.variables['coordx']
    y_1d = nc.variables['coordy']
    n_nodes = len(x_1d)
    ny = int(n_nodes/nx)
    print(f'nx = {nx}, ny = {ny}, n_nodes = {len(x_1d)}')
    assert(nx*ny==n_nodes)
    x = np.array(x_1d).reshape(ny,nx).T
    y = np.array(y_1d).reshape(ny,nx).T
    assert(np.abs(x[0,0]-x[0,1])<1e-8)

    # Extract the variables
    n_var = nc.variables['name_nod_var'].shape[0]
    # Decipher the order of the variables
    # for i in range(n_var):
    #     data_i = nc.variables['name_nod_var'][i].data
    #     filtered_data = [d.tobytes().decode() for d in data_i if d != b'--']
    #     print(''.join(filtered_data))
    y_d = nc.variables['vals_nod_var1'][-1,:].data.reshape(ny,nx).T
    P = nc.variables['vals_nod_var2'][-1,:].data.reshape(ny,nx).T
    omega = nc.variables['vals_nod_var3'][-1,:].data.reshape(ny,nx).T
    k = nc.variables['vals_nod_var4'][-1,:].data.reshape(ny,nx).T
    nu_t = nc.variables['vals_nod_var5'][-1,:].data.reshape(ny,nx).T
    U = nc.variables['vals_nod_var6'][-1,:].data.reshape(ny,nx).T
    V = nc.variables['vals_nod_var7'][-1,:].data.reshape(ny,nx).T

    Uref = 1.0
    rhoref = 1.0
    qref = 0.5*rhoref*Uref**2
    #Pref = 0.0
    Pref = P[3,0]
    # Re_L = 2e6
    # L = 1
    # nu = Uref*L/Re_L

    # Compute del_99 at x = -0.8
    x_target = -0.8
    idx_x = np.abs(x[:,0] - x_target).argmin()
    U_x = U[idx_x, :]
    V_x = V[idx_x, :]
    P_x = P[idx_x, :]
    Po = P_x + 0.5*rhoref*(U_x**2+V_x**2)
    Po_ref = np.max(Po)
    UI2_x = 2.0/rhoref*(Po_ref-P_x)-V_x**2
    indices_threshold = np.where(U_x**2 >= 0.99*UI2_x)[0]
    if len(indices_threshold) > 0:
        boundary_layer_thickness = y[idx_x,indices_threshold[0]]
    else:
        boundary_layer_thickness = 0.0
    print("GFM del99 at x = -0.8:", boundary_layer_thickness)
    # indices_threshold = np.where(U_x >= 0.99 * U_x[-1])[0]
    # if len(indices_threshold) > 0:
    #     boundary_layer_thickness = y[idx_x,indices_threshold[0]]
    # else:
    #     boundary_layer_thickness = 0.0
    # print("U_inf del99 at x = -0.8:", boundary_layer_thickness)
    # # Compute del_99 at x = -0.65
    x_target = -0.65
    idx_x = np.abs(x[:,0] - x_target).argmin()
    U_x = U[idx_x, :]
    V_x = V[idx_x, :]
    P_x = P[idx_x, :]
    Po = P_x + 0.5*rhoref*(U_x**2+V_x**2)
    Po_ref = np.max(Po)
    UI2_x = 2.0/rhoref*(Po_ref-P_x)-V_x**2
    indices_threshold = np.where(U_x**2 >= 0.99*UI2_x)[0]
    if len(indices_threshold) > 0:
        boundary_layer_thickness = y[idx_x,indices_threshold[0]]
    else:
        boundary_layer_thickness = 0.0
    print("GFM del99 at x = -0.65:", boundary_layer_thickness)
    # indices_threshold = np.where(U_x >= 0.99 * U_x[-1])[0]
    # if len(indices_threshold) > 0:
    #     boundary_layer_thickness = y[idx_x,indices_threshold[0]]
    # else:
    #     boundary_layer_thickness = 0.0
    # print("U_inf del99 at x = -0.65:", boundary_layer_thickness)

#fig.tight_layout()
#plt.show()

nc.close()
