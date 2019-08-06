import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from netCDF4 import Dataset


#---------------------------------------------------------------------------------------------------


# User provided input
# -------------------

# Date-time to plot
datetime = '20190401_0000z'

# Root path
root_path = '/gpfsm/dnb31/drholdaw/Victor/'

# Grid index to plot
i_index = 710
j_index = 2500

# Processor layout used for Jacobians
nx = 12
ny = nx*6

# Cube grid size
cube = 720


#---------------------------------------------------------------------------------------------------


# Paths
# -----

# Path to the Profiles and Jacobians training data
path_pro = root_path+'/TrainingData/'+datetime+'/Profiles/'
path_jac = root_path+'/TrainingData/'+datetime+'/Jacobians/'


# Checks
# ------

if i_index > cube or i_index < 1:
    print('i_index must be between 1 and '+str(cube))
    exit()
if j_index > 6*cube or j_index < 1:
    print('j_index must be between 1 and '+str(6*cube))
    exit()


# Domain decomposition for Jacobian calculation
# ---------------------------------------------

# Number of grid points per processor
ngp = 720/nx

# Determine is/js in file name for Jacobian
ifilename = i_index-(i_index-1)%ngp
jfilename = j_index-(j_index-1)%ngp

# Position within the file
ifile = i_index - ifilename + 1
jfile = j_index - jfilename + 1

# String verions
sisfile = str(ifilename).zfill(4)
sjsfile = str(jfilename).zfill(4)


# Open NetCDF file for reading
# ----------------------------

file_proin = path_pro+'/f522_dh.trainingdata_in.lcv.'+datetime+'.nc4'
file_prout = path_pro+'/f522_dh.trainingdata_out.lcv.'+datetime+'.nc4'
file_jacob = path_jac+'/jacobian_'+datetime+'_is'+sisfile+'_js'+sjsfile+'.nc4'

print('Files')
print(file_proin)
print(file_prout)
print(file_jacob)

fh_proin = Dataset(file_proin)
fh_prout = Dataset(file_prout)
fh_jacob = Dataset(file_jacob)


# Get lats and lons and check for equivalence
# -------------------------------------------

lon_pro = fh_proin.variables['lons'][j_index,i_index]
lon_jac = fh_jacob.variables['lons'][jfile,ifile]

lat_pro = fh_proin.variables['lats'][j_index,i_index]
lat_jac = fh_jacob.variables['lats'][jfile,ifile]

if abs(lat_pro-lat_jac) > 1.0e-9:
    print('Training and Jacobian latitudes dont seem to match, exit')
    exit()
if abs(lon_pro-lon_jac) > 1.0e-9:
    print('Training and Jacobian longitudes dont seem to match, exit')
    exit()


# Get profiles of temperature, pressure etc from training
# -------------------------------------------------------

t   = np.zeros([72])
pl  = np.zeros([72])
q   = np.zeros([72])
qi  = np.zeros([72])
ql  = np.zeros([72])
o3  = np.zeros([72])
flx = np.zeros([73])

t [:] = fh_proin.variables['t' ][:,:,j_index,i_index]
pl[:] = fh_proin.variables['pl'][:,:,j_index,i_index]
q [:] = fh_proin.variables['q' ][:,:,j_index,i_index]
qi[:] = fh_proin.variables['qi'][:,:,j_index,i_index]
ql[:] = fh_proin.variables['ql'][:,:,j_index,i_index]
o3[:] = fh_proin.variables['o3'][:,:,j_index,i_index]

flx[:] = fh_prout.variables['flxu'][:,:,j_index,i_index] + \
         fh_prout.variables['flxd'][:,:,j_index,i_index]


# Get Jacobians of temperature, pressure etc from training
# -------------------------------------------------------

dflxdt_jacob   = np.zeros([72,72])
dflxdpl_jacob  = np.zeros([72,72])
dflxdq_jacob   = np.zeros([72,72])
dflxdqi_jacob  = np.zeros([72,72])
dflxdql_jacob  = np.zeros([72,72])
dflxdo3_jacob  = np.zeros([72,72])

dflxdt_jacob[:,:]  = fh_jacob.variables['dflxdt'][:,:,jfile,ifile]
dflxdpl_jacob[:,:] = fh_jacob.variables['dflxdpl'][:,:,jfile,ifile]
dflxdq_jacob[:,:]  = fh_jacob.variables['dflxdq'][:,:,jfile,ifile]
dflxdqi_jacob[:,:] = fh_jacob.variables['dflxdqi'][:,:,jfile,ifile]
dflxdql_jacob[:,:] = fh_jacob.variables['dflxdql'][:,:,jfile,ifile]
dflxdo3_jacob[:,:] = fh_jacob.variables['dflxdo3'][:,:,jfile,ifile]


# Plot the training profiles for this grid point
# ----------------------------------------------

levs72 = np.arange(1,73)
levs73 = np.arange(1,74)

def add_profile_sp(ax,field,title,ylevs):

    ax.plot(field, ylevs)
    ax.set_ylim(1, np.max(ylevs))
    ax.invert_yaxis()
    ax.title.set_text(title)


fig1, axs1 = plt.subplots(2, 4, figsize=(25,15))

add_profile_sp(axs1[0,0], t,  't',  levs72)
add_profile_sp(axs1[0,1], pl, 'pl', levs72)
add_profile_sp(axs1[0,2], q,  'q',  levs72)
add_profile_sp(axs1[1,0], qi, 'qi', levs72)
add_profile_sp(axs1[1,1], ql, 'ql', levs72)
add_profile_sp(axs1[1,2], o3, 'o3', levs72)
add_profile_sp(axs1[0,3], flx, 'flx', levs73)

plt.savefig('profiles_'+datetime+'_is'+sisfile+'_js'+sjsfile+'.png')


# Contour plots of Jacobian
# --------------------------

def add_contour_sp(fig,ax,field,title):

    maxf = np.max(np.abs(field))
    if maxf == 0.0:
      maxf = 1.0
    incf = 2*maxf/21
    clevs = np.arange(-maxf,maxf+incf,incf)

    img = ax.contourf(levs72, levs72, np.transpose(field), clevs, cmap = 'PRGn')
    ax.set_aspect('equal', 'box')
    ax.invert_yaxis()
    ax.title.set_text(title)
    fig.colorbar(img,ax=ax)

fig2, axs2 = plt.subplots(2, 3, figsize=(25,15))

add_contour_sp(fig2, axs2[0,0], dflxdt_jacob,  'dflx/dt' )
add_contour_sp(fig2, axs2[0,1], dflxdpl_jacob, 'dflx/dpl')
add_contour_sp(fig2, axs2[0,2], dflxdq_jacob,  'dflx/dq' )
add_contour_sp(fig2, axs2[1,0], dflxdqi_jacob, 'dflx/dqi')
add_contour_sp(fig2, axs2[1,1], dflxdql_jacob, 'dflx/dql')
add_contour_sp(fig2, axs2[1,2], dflxdo3_jacob, 'dflx/do3')

plt.savefig('jacobian_'+datetime+'_is'+sisfile+'_js'+sjsfile+'.png')
