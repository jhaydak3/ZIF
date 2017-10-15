#!/usr/bin/env python 
import numpy as np
import sys
import readFiles
import density_f

def main():
	read_structure()
	read_trajectories()

def wrapCoordinates(xyz):
	#takes in a point
	for i in range(0,2):
		thisFrac = xyz[i]
		if thisFrac < 0:
			while thisFrac < 0:
				thisFrac = thisFrac + 1
		elif thisFrac >=1:
			while thisFrac >=1:
				thisFrac = thisFrac - 1
		xyz[i] = thisFrac
	return xyz
	
def frac2cart(xyz,cell_params):
	a = cell_params[0]; b = cell_params[1]; c = cell_params[2]
	alpha = np.deg2rad(cell_params[3]); beta = np.deg2rad(cell_params[4]); gamma = np.deg2rad(cell_params[5])
	omega = a*b*c * np.sqrt(1 - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2 + 2 * np.cos(alpha)*np.cos(beta)*np.cos(gamma))
	T = np.zeros((3,3))
	T[0,0] = a
	T[0,1] = b*np.cos(gamma)
	T[0,2] = c*np.cos(beta)
	T[1,1] = b*np.sin(gamma)
	T[1,2] = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
	T[2,2] = omega / a / b/ np.sin(gamma)
	nxyz = np.dot(T,xyz)
	return nxyz



def read_trajectories():
	global xx
	global yy
	global zz
	fh = open(traj_file,"r")
	zz = np.array(zz,order='F')
	xx = np.array(xx,order='F')
	yy = np.array(yy,order='F')
	numAtoms = int(fh.readline().strip())
	_ = fh.readline()
	_ = fh.readline()
	stayInWhileLoop = True
	cntr = 0.0
	while stayInWhileLoop:
		print cntr
		cntr = cntr+1.0
		for ndx1 in range(1,numAtoms):
			thisLine = fh.readline()
			if not thisLine:
				stayInWhileLoop = False
				break
			_, x, y, z = thisLine.split()
			x = float(x); y = float(y); z = float(z)
			xyz = readFiles.cart2frac(np.array([x,y,z]), cell_params)
			xyz = wrapCoordinates(xyz)
			xyz = frac2cart(xyz,cell_params)
			x = xyz[0]; y = xyz[1]; z = xyz[2]	
			if plane_flag == "xy":
				density_f.addcountxy(zz,x,y,xx,yy,r_cutoff)
			elif plane_flag == "xz":
				density_f.addcountxz(yy,x,z,xx,zz,r_cutoff)
			elif plane_flag == "yz":
				density_f.addcountyz(xx,y,z,yy,zz,r_cutoff)
			else:
				print("Plane flag not recognized")
		thisLine1 = fh.readline() #numatoms
		thisLine2 = fh.readline() #second line
		thisLine3 = fh.readline() #COM of structure
		if cntr % 500 == 0: #save every 500 loops
			if plane_flag == "xy":
				np.savetxt(traj_file[:-4] + "_densityxy.txt",zz,fmt="%d")
			elif plane_flag == "xz":
				np.savetxt(traj_file[:-4] + "_densityxz.txt",yy,fmt="%d")
			elif plane_flag == "yz":
				np.savetxt(traj_file[:-4] + "_densityyz.txt",xx,fmt="%d")
			else:
				print("Plane flag not recognized")
		if not thisLine1 or not thisLine2 or not thisLine3:
			stayInWhileLoop = False
			break
			
			
	if plane_flag == "xy":
		np.savetxt(traj_file[:-4] + "_densityxy.txt",zz,fmt="%d")
	elif plane_flag == "xz":
		np.savetxt(traj_file[:-4] + "_densityxz.txt",yy,fmt="%d")
	elif plane_flag == "yz":
		np.savetxt(traj_file[:-4] + "_densityyz.txt",xx,fmt="%d")
	else:
		print("Plane flag not recognized")
	fh.close()
	
	
def read_structure():
	global xx
	global yy
	global zz
	fh = open(struct_file,"r")
	numAtoms = int(fh.readline().strip())
	_ = fh.readline() #junk
	atomlist = []; xlist = []; ylist = []; zlist = []
	for i in range(0,numAtoms):
		thisAtom, thisx, thisy, thisz = fh.readline().split()
		atomlist.append(thisAtom)
		xlist.append(float(thisx))
		ylist.append(float(thisy))
		zlist.append(float(thisz))
	
	Znlistx = []; Znlisty = []; Znlistz = []
	for i in range(0,numAtoms):
		if atomlist[i] == "Zn":
			Znlistx.append(xlist[i])
			Znlisty.append(ylist[i])
			Znlistz.append(zlist[i])
	#write ZnCords to file
	ZnCoords = np.vstack((Znlistx,Znlisty,Znlistz))
	np.savetxt(struct_file[:-4] + "_ZnCoords.txt",ZnCoords,fmt="%f")
			
	if plane_flag == "xy":
		minx = min(Znlistx)
		maxx = max(Znlistx)
		miny = min(Znlisty)
		maxy = max(Znlisty)
		x = np.linspace(minx,maxx,ngrid)
		y = np.linspace(miny,maxy,ngrid)
		np.savetxt(struct_file[:-4] + "_xpoints.txt",x,fmt="%f")
		np.savetxt(struct_file[:-4] + "_ypoints.txt",y,fmt="%f")
		xx, yy = np.meshgrid(x,y)
		zz = 0*xx

	elif plane_flag == "xz":
		minx = min(Znlistx)
		maxx = max(Znlistx)
		minz = min(Znlistz)
		maxz = max(Znlistz)
		x = np.linspace(minx,maxx,ngrid)
		z = np.linspace(minz,maxz,ngrid)
		np.savetxt(struct_file[:-4] + "_xpoints.txt",x,fmt="%f")
		np.savetxt(struct_file[:-4] + "_zpoints.txt",z,fmt="%f")
		xx, zz = np.meshgrid(x,z)
		yy = 0*xx					
	
	elif plane_flag == "yz":
		miny = min(Znlisty)
		maxy = max(Znlisty)
		minz = min(Znlistz)
		maxz = max(Znlistz)
		y = np.linspace(miny,maxy,ngrid)
		z = np.linspace(minz,maxz,ngrid)
		np.savetxt(struct_file[:-4] + "_ypoints.txt",y,fmt="%f")
		np.savetxt(struct_file[:-4] + "_zpoints.txt",z,fmt="%f")
		yy, zz = np.meshgrid(y,z)
		xx = 0*yy	
	else:
		print("Plane flag not recognized!")

	fh.close()
	
	
	


if __name__ == "__main__":
	#INPUTS
	
	global cell_params
	cell_params = (44.9219, 44.9219, 31.89326, 90,90,120)
	global r_cutoff #distance to check for an atomlist
	r_cutoff = .5
	global ngrid
	ngrid = 250 #constructs an NxN grid to use
	global struct_file
	struct_file = "/nv/pj1/jhaydak3/Bash_Scripting/ZIF-7_DFTUNTYPED.xyz"
	global traj_file
	traj_file = "/nv/pj1/jhaydak3/FORCE_FIELD_MAKER/lmp_files/BzIm-ZIF-7-Hydrogen/MOLECULE_msd_traj.xyz"
	global plane_flag
	plane_flag = "xy" #xy, yz, xz
	global yy
	yy = None
	global zz
	zz = None
	global xx
	xx = None
	main()