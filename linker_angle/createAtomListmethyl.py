#!/usr/bin/env python 
import numpy as np


def createList():
	fh = open(xyz_file,"r")
	numAtoms = int(fh.readline())
	_ = fh.readline() #junk
	carbonList = []
	numCarbon = 0
	cartCords = np.empty([numAtoms,3])
	for i in range(0,numAtoms):
		thisLine = fh.readline()
		thisAtom, x, y, z = thisLine.split()
		cartCords[i,0] = x
		cartCords[i,1] = y
		cartCords[i,2] = z
		if thisAtom == "C":
			#print thisLine
			carbonList.append(i)
			numCarbon += 1
	carbonCords = np.empty([numCarbon, 3])
	#print carbonList
	for i in range(0,numCarbon):
		carbonCords[i,:] = cartCords[carbonList[i],:]
		#print carbonCords[i,:]
	fh.close()
	#carbonfracCords = cart2frac(carbonCords,params)
	#for i in range(0,numCarbon):
	#	print carbonfracCords[i,:]
	distMat = np.ones([numCarbon,numCarbon])*100
	for i in range(0,numCarbon):
		for j in range(i+1,numCarbon):
			dx = carbonCords[i,0] - carbonCords[j,0]
			dy = carbonCords[i,1] - carbonCords[j,1]
			dz = carbonCords[i,2] - carbonCords[j,2]
			thisdist = np.sqrt(dx**2 + dy**2 + dz**2)
			distMat[i,j] = thisdist
			distMat[j,i] = thisdist
	output_fh = open(output_string + '_atom_list.txt',"w")
	#read data file
	origin, A,B,C = readMinData()
	T = np.empty([3,3])
	T[:,0] = A; T[:,1] = B; T[:,2] = C
	invT = np.linalg.inv(T)
	print invT
	for i in range(0,numCarbon):
	#is it a prime carbon?
		minDist = np.min(distMat[i,:])
		minNdx = np.argmin(distMat[i,:])
		if minDist > 1.4: #then it is a primary carbon
			C2_1ndx = minNdx
			distMat[i,C2_1ndx] = distMat[i,C2_1ndx] + 100 #make this temporarily high so can call min again
			C2_2ndx = np.argmin(distMat[i,:])
			distMat[i,C2_1ndx] = distMat[i,C2_1ndx] - 100
			if distMat[i,C2_2ndx] < 3:
				distMat[i,C2_2ndx] = distMat[i,C2_2ndx] + 100
				distMat[i,C2_1ndx] = distMat[i,C2_1ndx] + 100
				C2_3ndx = np.argmin(distMat[i,:])
				distMat[i,C2_2ndx] = distMat[i,C2_1ndx] - 100
				distMat[i,C2_1ndx] = distMat[i,C2_1ndx] - 100
	

				#now check if all 3 within bounds
				#lowerBound = pad_factor
				#upperBound = 1-pad_factor
				#print i, C2_1ndx, C2_2ndx
				#print carbonfracCords[i,:]
				#print carbonfracCords[C2_1ndx,:]
				#print carbonfracCords[C2_2ndx,:]
				#checkifIn = np.all(carbonfracCords[i,:] < upperBound) and np.all(carbonfracCords[i,:] > lowerBound)
				#checkifIn = checkifIn and np.all(carbonfracCords[C2_1ndx,:] < upperBound) and np.all(carbonfracCords[C2_1ndx,:] > lowerBound)
				#checkifIn = checkifIn and np.all(carbonfracCords[C2_2ndx,:] < upperBound) and np.all(carbonfracCords[C2_2ndx,:] > lowerBound)
				#lowerBound = [limits[0]*pad_factor, limits[1]*pad_factor, limits[2]*pad_factor]
				#upperBound = [limits[0]*(1-pad_factor), limits[1]*(1-pad_factor), limits[2]*(1-pad_factor)]
				#checkifIn = np.all(carbonCords[i,:] > lowerBound) and np.all(carbonCords[i,:] < upperBound)
				#checkifIn = checkifIn and np.all(carbonCords[C2_1ndx,:] > lowerBound) and np.all(carbonCords[C2_1ndx,:] < upperBound)
				#checkifIn = checkifIn and np.all(carbonCords[C2_2ndx,:] > lowerBound) and np.all(carbonCords[C2_2ndx,:] < upperBound)
				thisPt = carbonCords[i,:] - origin
				thisC = np.dot(invT,thisPt)
				checkifIn = np.all(thisC < 1-pad_factor) and np.all(thisC > pad_factor)
				thisPt = carbonCords[C2_3ndx,:] - origin
				thisC = np.dot(invT,thisPt)
				checkifIn = checkifIn and np.all(thisC < 1-pad_factor) and np.all(thisC > pad_factor)
				thisPt = carbonCords[C2_2ndx,:] - origin
				thisC = np.dot(invT,thisPt)
				checkifIn = checkifIn and np.all(thisC < 1-pad_factor) and np.all(thisC > pad_factor)
				if checkifIn:
					output_fh.write(str(carbonList[i]) + ' ' + str(carbonList[C2_3ndx]) + ' ' + str(carbonList[C2_2ndx]) + '\n')
	output_fh.close()
				
			
			
			
			
		
def readMinData():
		fh = open(after_min_data_file,"r")
		for i in range(0,11):
			thisLine = fh.readline()
		xlo,xhi, _, _ = fh.readline().split()
		ylo,yhi, _, _ = fh.readline().split()
		zlo,zhi, _, _ = fh.readline().split()
		xy, xz, yz, _, _, _ = fh.readline().split()
		xlo = float(xlo); xhi = float(xhi)
		ylo = float(ylo); yhi = float(yhi)
		zlo = float(zlo); zhi = float(zhi)
		xy = float(xy); xz = float(xz); yz = float(xz)
		origin = np.array([xlo,ylo,zlo])
		A = np.array([xhi-xlo,0,0])
		B = np.array([xy, yhi-ylo, 0])
		C = np.array([xz,yz,zhi-zlo])
		return origin, A,B,C
	
		
		
		
		
	
def cart2frac(atomCord,cell_params):
	a = cell_params[0]; b = cell_params[1]; c = cell_params[2]
	alpha = np.deg2rad(cell_params[3]); beta = np.deg2rad(cell_params[4]); gamma = np.deg2rad(cell_params[5])
	omega = a*b*c * np.sqrt(1 - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2 + 2 * np.cos(alpha)*np.cos(beta)*np.cos(gamma))
	T = np.zeros((3,3))
	T[0,0] = 1/a
	T[0,1] = - np.cos(gamma) / a / np.sin(gamma)
	T[0,2] = b*c*(np.cos(alpha)*np.cos(gamma)-np.cos(beta))/ omega / np.sin(gamma)
	T[1,1] = 1 / b / np.sin(gamma)
	T[1,2] = a*c*(np.cos(beta)*np.cos(gamma)-np.cos(alpha)) / omega / np.sin(gamma)
	T[2,2] = a*b*np.sin(gamma) / omega
	atomCordFrac = np.dot(T,atomCord.T)
	atomCordFrac = atomCordFrac.T
	return atomCordFrac
	
def main():
	#INPUTS
	global xyz_file
	xyz_file = '/nv/pj1/jhaydak3/FORCE_FIELD_MAKER/lmp_files/ZIF-8/ZIF-8-Empty-300K/ZIF8_min.xyz'
	global after_min_data_file 
	after_min_data_file = '/nv/pj1/jhaydak3/FORCE_FIELD_MAKER/lmp_files/ZIF-8/ZIF-8-Empty-300K/AFTER_MIN_ROUTINE.data'
	global pad_factor
	pad_factor = .15
	global output_string
	output_string = 'ZIF8'
	#global params
	#params = (34.19463607, 34.22180789, 43.17457128,  82.37496683, 82.36807813, 109.47)
	#global limits
	#limits = (39.35821598, 29.93642006, 39.65451702)
	createList()
	
	
	
if __name__ == "__main__":
	main()