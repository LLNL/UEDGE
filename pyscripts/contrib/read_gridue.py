
def calcang(datarr):
	from numpy import pi,zeros,shape,arctan
# Find setups
	sets=list(datarr)
	setlen=len(sets)

# Loop through setups
	datarr["angpol"]=zeros(shape(datarr["rcc"]))
	datarr["angrad"]=zeros(shape(datarr["rcc"]))
	datarr["angpolraw"]=zeros(shape(datarr["rcc"]))
	datarr["angradraw"]=zeros(shape(datarr["rcc"]))
	for i in range(shape(datarr["rcc"])[0]):
		for j in range(shape(datarr["rcc"])[1]):
			lbx=0.5*(datarr["rll"][i,j]+datarr["rul"][i,j])	
			rbx=0.5*(datarr["rlr"][i,j]+datarr["rur"][i,j])	
			lbz=0.5*(datarr["zll"][i,j]+datarr["zul"][i,j])	
			rbz=0.5*(datarr["zlr"][i,j]+datarr["zur"][i,j])	
			dx=rbx-lbx
			dz=rbz-lbz
			if dx<0:
				sx=pi
			else:
				sx=0
			datarr["angpol"][i,j]=arctan(dz/dx)+sx
			datarr["angpolraw"][i,j]=arctan(dz/dx)

			ubx=0.5*(datarr["rul"][i,j]+datarr["rur"][i,j])
			lbx=0.5*(datarr["rll"][i,j]+datarr["rlr"][i,j])
			ubz=0.5*(datarr["zul"][i,j]+datarr["zur"][i,j])
			lbz=0.5*(datarr["zll"][i,j]+datarr["zlr"][i,j])
			dx=ubx-lbx
			dz=ubz-lbz
			if dx<0:
				sx=pi
			else:
				sx=0			
			datarr["angrad"][i,j]=arctan(dz/dx)+sx
			datarr["angradraw"][i,j]=arctan(dz/dx)





import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Open gridue as DF
df_gridue=pd.read_csv("../e2dscripts/gridue")

# First line contains geom. specifiers: extract and save them to variables
dims=list(df_gridue)[0].split()
dims=[int(x) for x in dims]
nx,ny,ixpt1,ixpt2,iysptrx=dims[0],dims[1],dims[2],dims[3],dims[4]
cells=(nx+2)*(ny+2)

# Convert to matrix
gridue=df_gridue.as_matrix() 
gridue=gridue[:-1,:] # Omit last line as it contains irrelevant data


# Data is read as strings. Split strings, replace decimal marker by E, convert tp
# float, and store in row vector
grid=[]
for i in range(len(gridue)):
	for j in range(3):
		grid.append(float(gridue[i][0].split()[j].replace("D","E")))

# Parse row vector in poloidal direction, starting from radially innermost position
data=[]
for i in range(len(grid)/(nx+2)):
	data.append(grid[i*(nx+2):(i+1)*(nx+2)])	
data=np.transpose(np.asarray(data))


# Each ny+2 rows make up a slice
slices=np.zeros((nx+2,ny+2,len(data.T)/(ny+2)))
for i in range(len(data.T)/(ny+2)):
	slices[:,:,i]=data[:,i*(ny+2):(i+1)*(ny+2)]

geo=dict()
pos=["cc","ll","lr","ul","ur"]
param=["r","z","psi_","br_","bz_","bpol_","bphi_","b_"]
for i in range(len(param)):
	for j in range(len(pos)):
		geo[param[i]+pos[j]]=slices[1:-1,1:-1,i*len(pos)+j]


grd=plt.figure()
ax1=grd.add_subplot(1,1,1)
point=1
depth=5
xs=[]
ys=[]
for j in range(5,10):
	for i in range(1,ny+1):
#	for x in range(1,nx):
#		plotx,ploty=[],[]
#		for j in [depth,depth+1,depth+2,depth+3,depth+4]:
#			plotx.append(grid[(point-1)*(cells*depth)+(nx+2)*i+1:(point-1)*(cells*depth)+(i+1)*(nx+2)-1][x])
#			ploty.append(grid[point*cells*depth+i*(nx+2)+1:point*depth*cells+(i+1)*(nx+2)-1][x])
#		plt.plot(plotx,ploty,"k-")

		
		xs.append((point-1)*(cells*j)+(nx+2)*i+1)
		xs.append((point-1)*(cells*j)+(i+1)*(nx+2)-1)
		ys.append(point*cells*j+i*(nx+2)+1)
		ys.append(point*j*cells+(i+1)*(nx+2)-1)
#		ax1.plot(grid[(point-1)*(cells*j)+(nx+2)*i+1:(point-1)*(cells*j)+(i+1)*(nx+2)-1],grid[point*cells*j+i*(nx+2)+1:point*j*cells+(i+1)*(nx+2)-1],"k-")
	#ax1.plot(grid[(point-1)*(cells*depth)+(nx+2)*i+1:(point-1)*(cells*depth)+(i+1)*(nx+2)-1],grid[point*cells*depth+i*(nx+2)+1:point*depth*cells+(i+1)*(nx+2)-1],"k-")



#	for x in range(1,nx):
#		plotx,ploty=[],[]
#		for j in [depth,depth+1,depth+2,depth+3,depth+4]:
#			plotx.append(grid[(point-1)*(cells*depth)+(nx+2)*i+1:(point-1)*(cells*depth)+(i+1)*(nx+2)-1][x])
#			ploty.append(grid[point*cells*depth+i*(nx+2)+1:point*depth*cells+(i+1)*(nx+2)-1][x])
#		plt.plot(plotx,ploty,"k-")

#	ax1.plot(grid[(point-1)*(cells*depth)+(nx+2)*i+1:(point-1)*(cells*depth)+(i+1)*(nx+2)-1],grid[point*cells*depth+i*(nx+2)+1:point*depth*cells+(i+1)*(nx+2)-1],"k-")


calcang(geo)



pol,rad=1,0
vix=geo["bpol_cc"]*pol*np.cos(geo["angpol"])+geo["bpol_cc"]*rad*np.cos(geo["angrad"])
viy=geo["bpol_cc"]*pol*np.sin(geo["angpol"])+geo["bpol_cc"]*rad*np.sin(geo["angrad"])
#vix=0.1*(pol*np.cos(geo["angpol"])+geo["bpol_cc"]*rad*np.cos(geo["angrad"]))
#viy=0.1*(pol*np.sin(geo["angpol"])+geo["bpol_cc"]*rad*np.sin(geo["angrad"]))




for i in range(np.shape(geo["rcc"])[0]):
	for j in range(np.shape(geo["rcc"])[1]):
		plotx,ploty=[],[]
		for k in ["ll","lr","ur","ul","ll"]:
			plotx.append(geo["r"+k][i,j])
			ploty.append(geo["z"+k][i,j])
#			plotx.append(slices[i,j,k])
#			ploty.append(slices[i,j,k+5])
		ax1.plot(plotx,ploty,"k-",linewidth=0.2)
#ax1.quiver(geo["rcc"],geo["zcc"],geo["br_cc"],0,scale=10)
#ax1.quiver(geo["rcc"],geo["zcc"],vix,viy,scale=10)
ax1.contourf(geo["rcc"],geo["zcc"],geo["bpol_cc"],np.linspace(0.015,0.73,20))
ax1.set_aspect("equal","datalim")
#ax1.plot(slices[1:-1,1:-1,0],slices[1:-1,1:-1,5],"k-")
#ax1.plot(slices[5,5,0:5],slices[5,5,5:10],"k-")
#ax1.plot(slices[5,5,0],slices[5,5,5],"ko")
#ax1.plot(slices[1:-1,-2,0],slices[1:-1,-2,5],"k-")
#ax1.plot(data[1:-1,ny],data[1:-1,5*(ny+2)+ny],"k-")
plt.suptitle("UE B_pol",fontsize=20)
grd.show()

"""

z=10
N=30
psi=plt.figure()
ax1=psi.add_subplot(1,1,1)

X=np.zeros((nx,ny*5))
Y=np.zeros((nx,ny*5))
Z=np.zeros((nx,ny*5))
for i in range(5):
	X[:,i*ny:(1+i)*ny]=slices[1:-1,1:-1,i]
	Y[:,i*ny:(1+i)*ny]=slices[1:-1,1:-1,i+5]
	Z[:,(i)*ny:(1+i)*ny]=slices[1:-1,1:-1,i+z]
CS=ax1.contourf(X,Y,Z,np.linspace(Z.min(),Z.max(),N))
cbar = psi.colorbar(CS)
cbar.ax.set_ylabel('Magnetic flux')
plt.title("Psi")
psi.show()


z=15
psi=plt.figure()
ax1=psi.add_subplot(1,1,1)

X=np.zeros((nx,ny*5))
Y=np.zeros((nx,ny*5))
Z=np.zeros((nx,ny*5))
for i in range(5):
	X[:,i*ny:(1+i)*ny]=slices[1:-1,1:-1,i]
	Y[:,i*ny:(1+i)*ny]=slices[1:-1,1:-1,i+5]
	Z[:,(i)*ny:(1+i)*ny]=slices[1:-1,1:-1,i+z]
CS=ax1.contourf(X,Y,Z,np.linspace(Z.min(),Z.max(),30))
cbar = psi.colorbar(CS)
cbar.ax.set_ylabel('Magnetic flux')
plt.title("B_r")
psi.show()


z=20
psi=plt.figure()
ax1=psi.add_subplot(1,1,1)

X=np.zeros((nx,ny*5))
Y=np.zeros((nx,ny*5))
Z=np.zeros((nx,ny*5))
for i in range(5):
	X[:,i*ny:(1+i)*ny]=slices[1:-1,1:-1,i]
	Y[:,i*ny:(1+i)*ny]=slices[1:-1,1:-1,i+5]
	Z[:,(i)*ny:(1+i)*ny]=slices[1:-1,1:-1,i+z]
CS=ax1.contourf(X,Y,Z,np.linspace(Z.min(),Z.max(),30))
cbar = psi.colorbar(CS)
cbar.ax.set_ylabel('Magnetic flux')
plt.title("B_z")
psi.show()


z=25
psi=plt.figure()
ax1=psi.add_subplot(1,1,1)

X=np.zeros((nx,ny*5))
Y=np.zeros((nx,ny*5))
Z=np.zeros((nx,ny*5))
for i in range(5):
	X[:,i*ny:(1+i)*ny]=slices[1:-1,1:-1,i]
	Y[:,i*ny:(1+i)*ny]=slices[1:-1,1:-1,i+5]
	Z[:,(i)*ny:(1+i)*ny]=slices[1:-1,1:-1,i+z]
CS=ax1.contourf(X,Y,Z,np.linspace(Z.min(),Z.max(),30))
cbar = psi.colorbar(CS)
cbar.ax.set_ylabel('Magnetic flux')
plt.title("B_pol")
psi.show()


z=30
psi=plt.figure()
ax1=psi.add_subplot(1,1,1)

X=np.zeros((nx,ny*5))
Y=np.zeros((nx,ny*5))
Z=np.zeros((nx,ny*5))
for i in range(5):
	X[:,i*ny:(1+i)*ny]=slices[1:-1,1:-1,i]
	Y[:,i*ny:(1+i)*ny]=slices[1:-1,1:-1,i+5]
	Z[:,(i)*ny:(1+i)*ny]=slices[1:-1,1:-1,i+z]
CS=ax1.contourf(X,Y,Z,np.linspace(Z.min(),Z.max(),30))
cbar = psi.colorbar(CS)
cbar.ax.set_ylabel('Magnetic flux')
plt.title("B_phi")
psi.show()



z=35
psi=plt.figure()
ax1=psi.add_subplot(1,1,1)

X=np.zeros((nx,ny*5))
Y=np.zeros((nx,ny*5))
Z=np.zeros((nx,ny*5))
for i in range(5):
	X[:,i*ny:(1+i)*ny]=slices[1:-1,1:-1,i]
	Y[:,i*ny:(1+i)*ny]=slices[1:-1,1:-1,i+5]
	Z[:,(i)*ny:(1+i)*ny]=slices[1:-1,1:-1,i+z]
CS=ax1.contourf(X,Y,Z,np.linspace(Z.min(),Z.max(),30))
cbar = psi.colorbar(CS)
cbar.ax.set_ylabel('Magnetic flux')
plt.title("B_tot")
psi.show()


# GRIDUE CRACKED
# Data is written in poloidal direction, from radially inner to outer position.
# Guard cell nodes are initialized to zeros - how about other values?
# First blocks are cell-center values, then SW,SE,NW,NE corners
# Parameters required for grid generation (in order) at each vertex
# -Radial vertex positions
# -Vertical vertex positions
# -Poloidal magnetic flux
# -Radial B-field strength
# -Vertical B-field strength
# -Poloidal B-field strength
# -Toroidal B-field strength
# -Total B-field strength

# Ordering:
# nx ny ixpt1(1) ixpt2(1) iysptrx1(1)
# Parameter-Vertex-Radial position-Poloidal position
#  Ordered into row vector, split at every third value




#gridue=np.asarray(grid).reshape(np.product(np.shape(gridue)))


#grid=np.zeros((cells,len(gridue)/cells))

#i=0
#while i <= 
#	for m < nx+2:
#		for n<ny+2:
#			grid[	


#plt.plot(grid[0:cells],grid[cells:2*cells])#plt.show()
#plt.show()
"""
