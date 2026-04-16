from __future__ import annotations

import vtk # for latex render
import sys
import os
import argparse 
import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
import pyvista

from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
elecDict = ParsedParameterFile('./.post/electricProperties')
# controlDict= ParsedParameterFile('./.post/controlDict')
 

img_folder = '.post'
# Ensure the directory exists
os.makedirs(img_folder, exist_ok=True)
# init reader
reader = pyvista.POpenFOAMReader("foam.foam")
mesh = reader.read()
# get some time space data
x_center = 0.5 * (mesh.bounds[0] + mesh.bounds[1])
x_length = np.abs(mesh.bounds[0] - mesh.bounds[1]) 
y_center =  0.5 * (mesh.bounds[2] + mesh.bounds[3])
y_length = np.abs(mesh.bounds[2] - mesh.bounds[3]) 
z_center = 0.5 * (mesh.bounds[4] + mesh.bounds[5])
z_length = np.abs(mesh.bounds[4] - mesh.bounds[5]) 
H = 2*y_length
times = reader.time_values 
 
 
parser = argparse.ArgumentParser(description="optional arguments")

### 1. Args for x slice plot 
parser.add_argument('--slicefield', type=str, default='Psi_E', help='Specified fieldto scope profile in z axis')
parser.add_argument('--slicexlocation', type=float, help='Specified x location to scope profile in z axis')
parser.add_argument('--slicetime', type=float, help='Specified time to scope profile in z axis')
args = parser.parse_args() 
 
# default handling
slicexlocation = x_center

# if user did not give a time to scope, pick the last time
if args.slicetime is not None:
    time = args.slicetime 
else:
    time = times[-1]


lambda_D = float(elecDict['lambda_D'])
print(lambda_D)

zeta = float(elecDict['PsiDict']['zeta'])
eps= float(elecDict['epsilon_0'])*float(elecDict['epsilon_r'])
mu = 1e-05*1000
 
plt.rcParams.update({
    "axes.labelsize": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14
})

def Psi_int_analytical(y):
    kappa = 1/lambda_D
    return zeta*np.cosh(kappa*y)/np.cosh(kappa*y_length)


def c_analytical(y,sign):
    k_bE = 1.38064852e-23;
    T = 300;
    e = 1.602176565e-19;
    return elecDict["c_bulk"] * (np.exp(  - sign * e * Psi_int_analytical(y) / (k_bE * T ) )   )


def U_analytical_free(y):
    kappa = 1/lambda_D
    E_x = float(elecDict['E_field_x'])
    Omega = eps * zeta / mu
    px =0
    print("px/dx =", px, "with L =", x_length) 
    return y_length*y_length*px/(2*mu)*(1-(y/y_length)*(y/y_length))\
            -(eps*E_x*zeta/mu )* (1- np.cosh(kappa*y)/np.cosh(kappa*y_length) )

def U_analytical_closed(y):
    kappa = 1/lambda_D
    E_x = float(elecDict['E_field_x'])
    Omega = eps * zeta / mu
    # return -(eps*E_x*zeta/mu )* (1- np.cosh(kappa*y)/np.cosh(kappa*y_length) )
    return -Omega * E_x * \
            ( -3/2*  \
                (1-(y/y_length)*(y/y_length)) \
                *(1- \
                    (np.tanh(kappa*y_length)/ \
                    (kappa*y_length))) \
                +(1-np.cosh(kappa*y)/np.cosh(kappa*y_length) ) \
            )
 


def main():
    # field =args.slicefield
    plot_Psi_E(times[-1],99)
    plot_concentration(times[-1],99,"cK")
    plot_concentration(times[-1],99,"cCl")
    plot_U(times[-1],99)
    plot_U_perp(times[-1],99)
    plot_p(times[-1],99)
 
    plot_U_free(times[-1],99)
 

    # plot_normalized_potential_int(50)
    # plot_c (times[-1])

# Extract a 2D slice (assuming the case is already 2D in one plane)
###### 0. Plot the 2D case 
def plot_2D_profile ():
    reader.set_active_time_value(time) ; mesh = reader.read()
    mesh_2D = mesh.slice(normal="y") 
    plotter = pyvista.Plotter(off_screen=True)
    plotter.add_mesh(mesh_2D, scalars=args.slicefield, cmap="magma", show_edges=True, lighting=False)  
    camera_position = (x_center, - y_center + 3* (z_center+ x_center), z_center)
    focal_point = (x_center, y_center, z_center)
    view_up = (0, 0, 1)
    plotter.camera_position = [camera_position, focal_point, view_up] 
    plotter.camera.roll = 0
    plotter.screenshot(os.path.join(img_folder, '2Dslice.png')) 
    plotter.show() 
    del mesh , mesh_2D



###### 1. plot a x-section slice of a field
def plot_Psi_E (time_in,L_percentage):
    #slice the mesh at time and space
    fig1, ax1 = plt.subplots() 
    indices = np.linspace(1, len(times) - 4, 5, dtype=int)
    plot_times = [times[i] for i in indices]    
    # plot_times =  times[1:]
    
    for t_val in plot_times:
        reader.set_active_time_value(t_val) ; mesh = reader.read()
        mesh_2D = mesh.slice(normal="z") 
        mesh_xslice = mesh_2D.slice(normal="x", origin=(x_center, 0, 0)) 

        fieldNameScalars =mesh_xslice["internalMesh"]["Psi_E"] # value
        y = mesh_xslice.get(0).cell_centers().points [:,1] # coords
        print (fieldNameScalars)
        print (Psi_int_analytical(y))
        # Plot scalar vs z
        ax1.plot(y/y_length,  fieldNameScalars/zeta, marker='o', label=f't={t_val:.2e}')
    ax1.plot( y/y_length,Psi_int_analytical(y)/zeta, marker='x',label='Analytical Solution')
    ax1.set_xlabel("$\\frac{z}{H/2}$")
    ax1.set_ylabel('$\\Psi_E/\\zeta$')
    ax1.set_title("")
    ax1.legend()
    ax1.grid(True)
    fig1.savefig(os.path.join(img_folder, '1Dxslice_Psi_E'+'@'+str(L_percentage)+'x.png'), dpi=200)  # Saves plot to file
    del mesh , mesh_2D, mesh_xslice


def plot_concentration (time_in,L_percentage,specieName):
    #slice the mesh at time and space
    fig1, ax1 = plt.subplots() 
    indices = np.linspace(1, len(times) - 2, 5, dtype=int)
    plot_times = [times[i] for i in indices]
    
    for t_val in plot_times:
        reader.set_active_time_value(t_val) ; mesh = reader.read()
        mesh_2D = mesh.slice(normal="z") 
        mesh_xslice = mesh_2D.slice(normal="x", origin=(x_center, 0, 0)) 

        fieldNameScalars =mesh_xslice["internalMesh"][specieName] # value
        y = mesh_xslice.get(0).cell_centers().points [:,1] # coords
        print (fieldNameScalars)
        # Plot scalar vs z
        ax1.plot(y/y_length,  fieldNameScalars/(elecDict["c_bulk"]), marker='o', label=f't={t_val:.2e}')
    if (specieName == "cK"):
        ax1.plot( y/y_length,c_analytical(y,+1)/(elecDict["c_bulk"]), marker='x',label='Analytical Solution')
    if (specieName =='cCl'):
        ax1.plot( y/y_length,c_analytical(y,-1)/(elecDict["c_bulk"]), marker='x',label='Analytical Solution')
    ax1.set_xlabel("$\\frac{z}{H/2}$")
    ax1.set_ylabel(rf"${specieName}/c_{{Bulk}}$")
    ax1.set_title("")
    ax1.legend()
    ax1.grid(True)
    fig1.savefig(os.path.join(img_folder, '1Dxslice_'+specieName+'@'+str(L_percentage)+'x.png'), dpi=200, bbox_inches='tight')  # Saves plot to file
    del mesh , mesh_2D, mesh_xslice

def plot_U (time_in,L_percentage):
    #slice the mesh at time and space
    fig1, ax1 = plt.subplots() 
    indices = np.linspace(1, len(times) - 4, 5, dtype=int)
    plot_times = [times[i] for i in indices]    
    # plot_times =  times[1:]
    Omega = eps * zeta / mu;     E_x = float(elecDict['E_field_x'])
    for t_val in plot_times:
        reader.set_active_time_value(t_val) ; mesh = reader.read()
        mesh_2D = mesh.slice(normal="z") 
        mesh_xslice = mesh_2D.slice(normal="x", origin=(x_center, 0, 0)) 

        print(mesh_xslice["internalMesh"])
        fieldNameScalars =mesh_xslice["internalMesh"]["U"] # value
        fieldNameScalars = fieldNameScalars [:, 0]
        y = mesh_xslice.get(0).cell_centers().points [:,1] # coords
        print (fieldNameScalars)
        # Plot scalar vs z
        ax1.plot(y/y_length,  fieldNameScalars/(Omega*E_x), marker='o', label=f't={t_val:.2e}')
    ax1.plot( y/y_length,U_analytical_closed(y)/(Omega*E_x), marker='x',label='Analytical Solution, closed channel')
    ax1.set_xlabel("$\\frac{z}{H/2}$")
    ax1.set_ylabel("$\\frac{U_x}{\\Omega*E_x}$")
    ax1.set_title("")
    ax1.legend()
    ax1.grid(True)
    fig1.savefig(os.path.join(img_folder, '1Dxslice_U_X'+'@'+str(L_percentage)+'x.png'), dpi=200, bbox_inches='tight')  # Saves plot to file
    del mesh , mesh_2D, mesh_xslice

def plot_U_perp (time_in,L_percentage):
    #slice the mesh at time and space
    fig1, ax1 = plt.subplots() 
    indices = np.linspace(1, len(times) - 4, 5, dtype=int)
    plot_times = [times[i] for i in indices]    
    # plot_times =  times[1:]
    Omega = eps * zeta / mu;     E_x = float(elecDict['E_field_x'])
    for t_val in plot_times:
        reader.set_active_time_value(t_val) ; mesh = reader.read()
        mesh_2D = mesh.slice(normal="z") 
        mesh_xslice = mesh_2D.slice(normal="x", origin=(x_center, 0, 0)) 

        print(mesh_xslice["internalMesh"])
        fieldNameScalars =mesh_xslice["internalMesh"]["U"] # value
        fieldNameScalars = fieldNameScalars [:, 1]
        y = mesh_xslice.get(0).cell_centers().points [:,1] # coords
        print (fieldNameScalars)
        # Plot scalar vs z
        ax1.plot(y/y_length,  fieldNameScalars/(Omega*E_x), marker='o', label=f't={t_val:.2e}')
    ax1.set_xlabel("$\\frac{z}{H/2}$")
    ax1.set_ylabel("$\\frac{U_y}{\\Omega*E_x}$")
    ax1.set_title("")
    ax1.legend()
    ax1.grid(True)
    fig1.savefig(os.path.join(img_folder, '1Dxslice_U_Y'+'@'+str(L_percentage)+'x.png'), dpi=200, bbox_inches='tight')  # Saves plot to file
    del mesh , mesh_2D, mesh_xslice

def plot_p(time_in,L_percentage):
    #slice the mesh at time and space
    fig1, ax1 = plt.subplots() 
    indices = np.linspace(1, len(times) - 4, 5, dtype=int)
    plot_times = [times[i] for i in indices]    
    # plot_times =  times[1:]
    Omega = eps * zeta / mu;     E_x = float(elecDict['E_field_x'])
    for t_val in plot_times:
        reader.set_active_time_value(t_val) ; mesh = reader.read()
        mesh_2D = mesh.slice(normal="z") 
        mesh_xslice = mesh_2D.slice(normal="x", origin=(x_center, 0, 0)) 

        print(mesh_xslice["internalMesh"])
        fieldNameScalars =mesh_xslice["internalMesh"]["p"] # value
        y = mesh_xslice.get(0).cell_centers().points [:,1] # coords
        print (fieldNameScalars)
        # Plot scalar vs z
        ax1.plot(y/y_length,  fieldNameScalars, marker='o', label=f't={t_val:.2e}')
    ax1.set_xlabel("$\\frac{z}{H/2}$")
    ax1.set_ylabel("$p$")
    ax1.set_title("")
    ax1.legend()
    ax1.grid(True)
    fig1.savefig(os.path.join(img_folder, '1Dxslice_p'+'@'+str(L_percentage)+'x.png'), dpi=200, bbox_inches='tight')  # Saves plot to file
    del mesh , mesh_2D, mesh_xslice

def plot_U_free (time_in,L_percentage):
    #slice the mesh at time and space
    fig1, ax1 = plt.subplots() 
    indices = np.linspace(1, len(times) - 4, 5, dtype=int)
    plot_times = [times[i] for i in indices]    
    # plot_times =  times[1:]
    Omega = eps * zeta / mu;     E_x = float(elecDict['E_field_x'])
    for t_val in plot_times:
        reader.set_active_time_value(t_val) ; mesh = reader.read()
        mesh_2D = mesh.slice(normal="z") 
        mesh_xslice = mesh_2D.slice(normal="x", origin=(x_center, 0, 0)) 

        print(mesh_xslice["internalMesh"])
        fieldNameScalars =mesh_xslice["internalMesh"]["U"] # value
        fieldNameScalars = fieldNameScalars [:, 0]
        y = mesh_xslice.get(0).cell_centers().points [:,1] # coords
        print (fieldNameScalars)
        # Plot scalar vs z
        ax1.plot(y/y_length,  fieldNameScalars/(Omega*E_x), marker='o', label=f't={t_val:.2e}')
    ax1.plot( y/y_length,U_analytical_free(y)/(Omega*E_x), marker='x',label='Analytical Solution, free channel')
    ax1.set_xlabel("$\\frac{z}{H/2}$")
    ax1.set_ylabel("$\\frac{U_x}{\\Omega*E_x}$")
    ax1.set_title("")
    ax1.legend()
    ax1.grid(True)
    fig1.savefig(os.path.join(img_folder, '1Dxslice_U_X_free'+'@'+str(L_percentage)+'x.png'), dpi=200, bbox_inches='tight')  # Saves plot to file
    del mesh , mesh_2D, mesh_xslice
 



if __name__ == "__main__":
    main()
