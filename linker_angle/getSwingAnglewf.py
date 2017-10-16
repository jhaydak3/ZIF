#!/usr/bin/env python 
import numpy as np
import swingAngle

 
def readDihedralsFromxyz(cordList):
	#fh = open(xyz_file,"r")
	#first read in reference atom positions ( of minimized structure)
	ref_fh = open(minimized_xyz_file,"r")
	numAtoms = getNumAtoms(ref_fh)
	ref_fh.close()
	#set up array for storing angle vs time for each linker
	numCord = len(cordList)
	timePlot = np.empty([numStepsToRead, numCord+1])
	ndxArray = np.array(cordList,order='F')
	numLinkers = ndxArray.shape[0]
	numAngles = numStepsToRead*numLinkers
	thetaList = np.empty([numAngles])
	thetaList = np.array(thetaList,order='F')
	timePlot = np.array(timePlot,order='F')
	swingAngle.makeanglelist(timePlot,xyz_file,minimized_xyz_file,numAtoms,thetaList,ndxArray)
	
	#for i in range(0,numStepsToRead):
		
		# thisDict, timestep = makeCordDict(fh,ndxToGrabList)
		# timePlot[i,0] = timestep
		# j = 1
		# if i % 500 == 0:
			# print i
		# for thisList in cordList:
			# C1 = thisDict[thisList[0]]
			# C2 = thisDict[thisList[1]]
			# C3 = thisDict[thisList[2]]
			# C = np.array(C2) + np.array(C3)
			# C = C/2
			# v1 = C - C1;
			
			# C1ref = refDict[thisList[0]]
			# C2ref = refDict[thisList[1]]
			# C3ref = refDict[thisList[2]]
			# Cref = np.array(C2ref) + np.array(C3ref)
			# Cref = Cref/2
			# v1ref = Cref - C1ref;
			
			# thisTheta = findAnglebtwnVecs(v1,v1ref)
			# thisTheta = np.rad2deg(thisTheta)
			# thetaList.append(thisTheta)
			# timePlot[i,j] = thisTheta
			# j = j + 1


	#fh.close()

	np.savetxt(save_file_name + "time_series.txt",timePlot,fmt="%f")
	return thetaList
	


	
def getNumAtoms(fh):
	numAtoms = fh.readline().strip()
	numAtoms = int(numAtoms)
	return numAtoms
		
		
	
	
	
def getList():
	#creates a list of dihedral coords of the form C1 C2 C2
	fh = open(atom_list,"r")
	thisLine = fh.readline()
	cordList = [];
	while(thisLine):	
		C1, C2, C3 = thisLine.split()
		C1 = int(C1); C2 = int(C2); C3 = int(C3)
		thisList = [C1, C2, C3]
		cordList.append(thisList)
		thisLine = fh.readline()
	return cordList
	fh.close()
		
def makeBins(thetaList,bin_size,save_file_name):
	lower_theta = -180
	upper_theta = 180
	bins = np.arange(lower_theta,upper_theta,bin_size)
	thisDict = {}
	total = 0.0
	ndx = 0
	#print thetaList
	for theta in thetaList:
	#	theta = np.abs(theta)
	#	if theta <0:
	#		theta = theta + 180
	#	if theta >=90:
	#		theta = 180 - theta
		theta = np.nan_to_num(theta)
		thetaList[ndx] = theta
		ndx = ndx + 1
	for bin in bins:
		thisDict[bin] = 0
	for theta in thetaList:
		remainder = theta % bin_size
		theta = theta - remainder
		thisDict[theta] += 1
		total += 1
	#now we have an unsorted dict
	sortedKeys = sorted(thisDict.keys())
	binList = []
	k = 8.3144598 # J/mol K
	k = k/1000 #kJ/mol K
	for key in sortedKeys:
		thisCount = thisDict[key]
		thisProb= thisCount/total	
		if thisProb > 0:
			thisF = -k*TEMP*np.log(thisProb)
		else:
			thisF = 0
		binList.append([key,thisProb,thisF])
	binList = np.array(binList)
	np.savetxt(save_file_name + "_angle_results.txt",binList,fmt="%f")
		

def main():
	#INPUTS
	global xyz_file
	xyz_file = '/nv/pj1/jhaydak3/FORCE_FIELD_MAKER/lmp_files/neb_empty_300/MOVIE_PROD.xyz'
	global minimized_xyz_file #used for reference angles
	minimized_xyz_file = '/nv/pj1/jhaydak3/FORCE_FIELD_MAKER/lmp_files/neb_empty_300/neb_empty_min.xyz'
	global atom_list
	atom_list = '/nv/pj1/jhaydak3/FORCE_FIELD_MAKER/lmp_files/neb_empty_300/NEB_atom_list.txt'
	global numStepsToRead
	numStepsToRead = 15000
	global bin_size
	bin_size = .25 #degrees
	global save_file_name 
	save_file_name= 'NEbtesting'
	global TEMP
	TEMP = 308.15 #Kelvin
	cordList = getList()
	thetaList = readDihedralsFromxyz(cordList)
	makeBins(thetaList,bin_size,save_file_name)

if __name__ == "__main__":
	main()