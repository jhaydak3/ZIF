import numpy as np
import pdb
import os.path

def readLatParamsFromxyz(fh):
	#reads lattice parameters from xyz file if they are in the second line
	rfile = open(fh,"r")
	thisLine = rfile.readline()
	thisLine = rfile.readline()
	a, b, c, alpha, beta, gamma = thisLine.split()
	rfile.close()
	return a,b,c,alpha,beta,gamma
	

def readxyz(fh):
#Reads an xyz file, returns a list of atoms and matrix of xyz coordinates
	rfile = open(fh,"r")
	numAtoms = int(rfile.readline())
	rfile.readline() # read the blank line
	atomList = []
	atomCord = np.empty([numAtoms,3])
	for i in range(0,numAtoms):
		thisAtom,thisx,thisy,thisz = rfile.readline().split()
		atomList.append(thisAtom)
		atomCord[i,0] = float(thisx)
		atomCord[i,1] = float(thisy)
		atomCord[i,2] = float(thisz)
	rfile.close()
	return atomList, atomCord
	
def readcifFile(fh):
# reads a .cif file, returns list of atoms and xyz coordinates
	sfile = open(fh,"r")
	foundAll = False; foundA = False; foundB = False; foundC = False
	foundAlpha = False; foundBeta = False; foundGamma = False
	atomCord = np.empty([0,3])
	#Read until get to start of coordinates
	while not foundAll:
		thisLine = sfile.readline()
		if "_cell_length_a" in thisLine:
			_junk, cella = thisLine.split()
			cella = float(cella)
			foundA = True
		elif "_cell_length_b" in thisLine:
			_junk, cellb = thisLine.split()
			cellb = float(cellb)
			foundB = True
		elif "_cell_length_c" in thisLine:
			_junk, cellc = thisLine.split()
			cellc = float(cellc)
			foundC = True
		elif "_cell_angle_alpha" in thisLine:
			_junk, alpha = thisLine.split()
			alpha = float(alpha)
			foundAlpha = True
		elif "_cell_angle_beta" in thisLine:
			_junk, beta = thisLine.split()
			beta = float(beta)
			foundBeta = True
		elif "_cell_angle_gamma" in thisLine:
			_junk, gamma = thisLine.split()
			gamma = float(gamma)
			foundGamma = True
		foundAll = foundA and foundB and foundC and foundAlpha and foundBeta and foundGamma
	return cella, cellb, cellc, alpha, beta, gamma
	
def writexyzFile(atomCord,atomListU,fh,label = "cheese"):
	wfile = open(fh,"w")
	wfile.write(str(atomCord.shape[0]) + "\n")
	if label == "cheese":
		wfile.write("No unit cell parameter/angle info available." + "\n")
	else:
		wfile.write(str(label[0]) + " " + str(label[1]) + " " + str(label[2]))
		wfile.write(" " + str(label[3]) + " " + str(label[4]) + " " + str(label[5]) + "\n")
	for i in range(0,len(atomListU)):
		wfile.write(atomListU[i] + "\t")
		wfile.write(str(atomCord[i,0]) + "\t")
		wfile.write(str(atomCord[i,1]) + "\t")
		wfile.write(str(atomCord[i,2]) + "\n")
	wfile.close()
	return 0
	
def replaceSiwZn(atomList):
	natomList = atomList
	numAtoms = len(natomList)
	for i in range(0,numAtoms):
		if natomList[i] == "Si":
			natomList[i] = "Zn"
	return natomList
	
def expandStructure(atomCord,sf):
	T = np.array([[sf, 0, 0],[0, sf ,0],[0, 0, sf]])
	natomCord = np.dot(T,atomCord.T)
	natomCord = natomCord.T
	return natomCord

def calcMinZnZnDist(atomCord,atomList):
# finds minimum Zn-Zn bond distance for scaling
	ZnPos = []
	for i in range(0,len(atomList)):
		if atomList[i] == "Zn":
			ZnPos.append(i)
	numZnatoms = len(ZnPos)
	distZnZn = 100*np.ones([numZnatoms,numZnatoms])
	for i in range(0,numZnatoms):
		for j in range(i+1,numZnatoms):
			dx = atomCord[ZnPos[i],0] - atomCord[ZnPos[j],0]
			dy = atomCord[ZnPos[i],1] - atomCord[ZnPos[j],1]
			dz = atomCord[ZnPos[i],2] - atomCord[ZnPos[j],2]
			distZnZn[i,j] = np.sqrt(dx**2+dy**2+dz**2)
	#np.savetxt("ZnZndist.txt",distZnZn, fmt ="%2.3f")
	#need to find min, but 0
	minDist = np.amin(distZnZn)
	return minDist
	
def putLinkerIn(atomListS,atomListL,atomCordS,atomCordL):
#replaces O atoms with linker
	x,y,z = calcLinkerCenter(atomCordL,atomListL)
	#nLinkerCord = translateCord(atomCordL,-x,-y,-z) #linker centered at 0
	nLinkerCord = noMoreNaughtyLinker(atomCordL,atomListL)
	#writexyzFile(nLinkerCord,atomListL,"linker_zero.xyz")
	#first count number of atoms in things
	numSAtoms = len(atomListS)
	numLAtoms  = len(atomListL)
	#count # of O atoms, doing this so we know how big to make a matrix of the new 
	# coordinates. each O atom is removed and replaced with the Linker atoms
	numOatoms = 0
	for i in range(0,numSAtoms):
		if atomListS[i] == "O":
			numOatoms = numOatoms + 1
	newnumSAtoms = numSAtoms - numOatoms + numOatoms*numLAtoms
	newatomCordS = np.empty([newnumSAtoms,3])
	#go through the atoms again, this time building new coordinates
	newatomListS = []
	j = 0 #index to loop through newatomCordS
	for i in range(0,numSAtoms):
		if atomListS[i] != "O": #not an oxygen, add it in
			newatomListS.append(atomListS[i])
			newatomCordS[j,0] = atomCordS[i,0]
			newatomCordS[j,1] = atomCordS[i,1]
			newatomCordS[j,2] = atomCordS[i,2]
			j = j+1
		else: #add in a linker at the coordinates of the O atom
			x = atomCordS[i,0]
			y = atomCordS[i,1]
			z = atomCordS[i,2]
			thisLinkerCord = translateCord(nLinkerCord,x,y,z)
			#get positions of the two closest Zn atoms
			Zn1,Zn2 = findClosestZn(x,y,z,atomCordS,atomListS)
			#check to see if the closest Zn atom is well behaved - 
			# right now i don't want to deal with atoms on the edge 
			isPlaceholder1 = np.all(Zn1 == 100)
			isPlaceholder2 = np.all(Zn2 == 100)
			#right now fuck the boundary atoms
			if not isPlaceholder1 and not isPlaceholder2:
				O_pos = np.array([x,y,z])
				thisLinkerCord = rotateThisBitch(Zn1,Zn2,O_pos,thisLinkerCord,atomListL)
			#now loop through the linker atoms and add them in to the new structure list and coordinates
			for k in range(0,numLAtoms):
				newatomListS.append(atomListL[k])
				newatomCordS[j,0] = thisLinkerCord[k,0]
				newatomCordS[j,1] = thisLinkerCord[k,1]
				newatomCordS[j,2] = thisLinkerCord[k,2]
				j = j + 1
	#writexyzFile(newatomCordS,newatomListS,"test_linker_in_not_rotated.xyz")		
	#writexyzFile(newatomCordS,newatomListS,"test_linker_rotated.xyz")		
	return newatomCordS,newatomListS

def noMoreNaughtyLinker(thisLinkerCord,atomListL):
#takes in the linker coordinates and translates/rotates it so that it is lying on the xy plane with its center at the origin
#start out by finding the C and N indices
	N_ndx = []; C_ndx = []
#ASSUMES NO CARBON GROUPS - IF THERE ARE THIS MAY NOT WORK
#WARNING
	numAtoms = len(atomListL)
	for i in range(0,numAtoms):
		if atomListL[i] == "C":
			C_ndx.append(i)
		elif atomListL[i] == "N":
			N_ndx.append(i)
	#translate the linker so that one of the Nitrogens is on the origin
	Nx = thisLinkerCord[N_ndx[0],0]; Ny = thisLinkerCord[N_ndx[0],1]; Nz = thisLinkerCord[N_ndx[0],2]
	thisLinkerCord = translateCord(thisLinkerCord,-Nx,-Ny,-Nz)
	#rotate second N about Z axis so that it lies on xz plane
	n1 = np.array([0,1,0])
	v1 = thisLinkerCord[N_ndx[1],:]
	v2 = np.array([thisLinkerCord[N_ndx[1],0], thisLinkerCord[N_ndx[1],1], 0])
	n2 = np.cross(v1,v2)
	theta = findAnglebtwnVecs(n1,n2)
	thisLinkerCord = rotateAboutZAxis(thisLinkerCord,theta)
	#make sure rotation in the right direction
	y_check = thisLinkerCord[N_ndx[1],1]
	if np.abs(y_check) > 1e-5:
		thisLinkerCord = rotateAboutZAxis(thisLinkerCord,-2*theta)
	#rotate second N about x axis so that it lies on the x axis
	n1 = np.array([0,0,1])
	v1 = thisLinkerCord[N_ndx[1],:]
	v2 = np.array([1,0,thisLinkerCord[N_ndx[1],2]])
	n2 = np.cross(v1,v2)
	theta = findAnglebtwnVecs(n1,n2)
	thisLinkerCord = rotateAboutXAxis(thisLinkerCord,theta)
	#make sure rotation in the right direction
	z_check = thisLinkerCord[N_ndx[1],2]
	if np.abs(z_check) > 1e-5:
		thisLinkerCord = rotateAboutXAxis(thisLinkerCord,-2*theta)
	#finally, get entire molecule to be on xy plane - don't need to use normals, can consider
	#2-d projection for calculating last angle
	v1 = np.array([0,1,0])
	v2 = np.array([0, thisLinkerCord[C_ndx[0],1], thisLinkerCord[C_ndx[0],2]])
	theta = findAnglebtwnVecs(v1,v2)
	thisLinkerCord = rotateAboutXAxis(thisLinkerCord,theta)
	z_check = thisLinkerCord[C_ndx[0],2]
	if np.abs(z_check) > 1e-5:
		thisLinkerCord = rotateAboutXAxis(thisLinkerCord,-2*theta)
	#move it so that the center is in the origin
	cx,cy,cz = calcLinkerCenter(thisLinkerCord,atomListL)
	thisLinkerCord = translateCord(thisLinkerCord,-cx,-cy,-cz)
	return thisLinkerCord
	
def rotateThisBitch(Zn1,Zn2,O_pos,thisLinkerCord,atomListL):
	N_ndx = []; C_ndx = []
	#ASSUMES NO CARBON GROUPS - IF THERE ARE THIS MAY NOT WORK
	#WARNING
	numAtoms = len(atomListL)
	for i in range(0,numAtoms):
		if atomListL[i] == "C":
			C_ndx.append(i)
		elif atomListL[i] == "N":
			N_ndx.append(i)
	#first, make it so the O atom is at the origin
	ZnCord = np.append([Zn1],[Zn2],axis=0)
	Ox = O_pos[0]; Oy = O_pos[1]; Oz = O_pos[2]
	ZnCord = translateCord(ZnCord,-Ox,-Oy,-Oz)
	thisLinkerCord = translateCord(thisLinkerCord,-Ox,-Oy,-Oz)
	#want to get the Zn's to be on xy plane
	#rotate Zn1 about x axis so that it is on xy plane
	#pdb.set_trace()
	v1 = np.array([0,1,0])
	v2 = np.array([0,ZnCord[0,1],ZnCord[0,2]])
	theta_xaxis = findAnglebtwnVecs(v1,v2)
	ZnCord = rotateAboutXAxis(ZnCord,theta_xaxis)
	z_check = ZnCord[0,2]
	if np.abs(z_check) > 1e-4:
		theta_xaxis = -theta_xaxis
		ZnCord = rotateAboutXAxis(ZnCord,2*theta_xaxis)
	#rotate Zn1 about z axis so that it is on x axis
	v1 = np.array([1,0,0])
	v2 = ZnCord[0,:]
	theta_zaxis = findAnglebtwnVecs(v1,v2)
	ZnCord = rotateAboutZAxis(ZnCord,theta_zaxis)
	y_check = ZnCord[0,1]
	if np.abs(y_check) > 1e-4:
		theta_zaxis = -theta_zaxis
		ZnCord = rotateAboutZAxis(ZnCord,2*theta_zaxis)
	#rotate Zn2 about x axis so that it is on the xy plane
	theta_xaxis2 = 0
	if np.abs(ZnCord[1,2]) > 1e-4:
		#pdb.set_trace()
		n1 = np.array([0,0,1])
		v1 = ZnCord[1,:]
		#v2 = np.array([ZnCord[1,0],ZnCord[1,1],0])
		v2 = np.array([0,ZnCord[1,1],ZnCord[1,2]])
		n2 = np.cross(v1,v2)
		theta_xaxis2 = findAnglebtwnVecs(n1,n2)
		ZnCord = rotateAboutXAxis(ZnCord,theta_xaxis2)
		z_check = ZnCord[1,2]
		if np.abs(z_check) > 1e-4:
			theta_xaxis2 = -theta_xaxis2
			ZnCord = rotateAboutXAxis(ZnCord,2*theta_xaxis2)
	
	#identify the primary carbon
	N1 = thisLinkerCord[N_ndx[0],:]; N2 = thisLinkerCord[N_ndx[1],:]
	C_prim_ndx = "hi" #error checking later on
	for i in range(0,len(C_ndx)):
		#assume that C-N bond is < 1.5 Ang
		thisC = thisLinkerCord[C_ndx[i],:]
		dx = thisC[0] - N1[0]; dy = thisC[1] - N1[1]; dz = thisC[2] - N1[2]
		N1dist = np.sqrt(dx**2+dy**2+dz**2)
		dx = thisC[0] - N2[0]; dy = thisC[1] - N2[1]; dz = thisC[2] - N2[2]
		N2dist = np.sqrt(dx**2+dy**2+dz**2)
		#print N1dist,N2dist
		if N1dist < 1.5 and N2dist < 1.5:
			C_prim_ndx = C_ndx[i]
			break
	if C_prim_ndx == "hi":
	#with a good linker this should never happen
		raise ValueError('Primary Carbon not found!')
		
	#rotate the linker so that it N atoms lined up with Zn atoms
	#pdb.set_trace()
	v1 = thisLinkerCord[C_prim_ndx,:]
	v2 = (ZnCord[0,:] + ZnCord[1,:]) / 2 # midpoint of Zn atoms
	theta_zaxis2= findAnglebtwnVecs(v1,v2)
	thisLinkerCord = rotateAboutZAxis(thisLinkerCord,theta_zaxis2)
	#check to make sure in right direction
	v1 = thisLinkerCord[C_prim_ndx,:]
	theta_check = findAnglebtwnVecs(v1,v2)
	if theta_check > 1e-3:
		theta_zaxis2 = -theta_zaxis2
		thisLinkerCord = rotateAboutZAxis(thisLinkerCord,2*theta_zaxis2)
	
	testS = np.concatenate((thisLinkerCord,ZnCord),axis = 0)
	atomListNew = atomListL[:]
	atomListNew.append('Zn')
	atomListNew.append('Zn')
	# fh = "test_frag"
	# thisfh = fh
	# if os.path.exists(fh + ".xyz"):
		# k = 1
		# while os.path.exists(thisfh+".xyz"):
			# thisfh = fh + str(k)
			# k = k + 1
	# writexyzFile(testS,atomListNew,thisfh+'.xyz')
	#now undo everything
	thisLinkerCord = rotateAboutXAxis(thisLinkerCord,-theta_xaxis2)
	thisLinkerCord = rotateAboutZAxis(thisLinkerCord,-theta_zaxis)
	thisLinkerCord = rotateAboutXAxis(thisLinkerCord,-theta_xaxis)
	thisLinkerCord = translateCord(thisLinkerCord,Ox,Oy,Oz)
	ZnCord = rotateAboutXAxis(ZnCord,-theta_xaxis2)
	ZnCord = rotateAboutZAxis(ZnCord,-theta_zaxis)
	ZnCord = rotateAboutXAxis(ZnCord,-theta_xaxis)
	ZnCord = translateCord(ZnCord,Ox,Oy,Oz)
	#testS = np.concatenate((thisLinkerCord,ZnCord),axis = 0)
	#atomListNew = atomListL[:]
	#atomListNew.append('Zn')
	#atomListNew.append('Zn')
	#fh = "test_frag"
	#thisfh = fh
	# if os.path.exists(fh + ".xyz"):
		# k = 1
		# while os.path.exists(thisfh+".xyz"):
			# thisfh = fh + str(k)
			# k = k + 1
	#writexyzFile(testS,atomListNew,thisfh+'.xyz')
	
	return thisLinkerCord

	
def rotateAboutZAxis(atomCord, theta):
	atomCord = atomCord.T
	Rz = np.array([[np.cos(theta), -np.sin(theta), 0], [ np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
	atomCord = np.dot(Rz,atomCord)
	return atomCord.T

def rotateAboutYAxis(atomCord,theta):
	atomCord = atomCord.T
	Ry = np.array([[np.cos(theta), 0 , np.sin(theta)], [ 0, 1, 0], [-np.sin(theta), 0, np.cos(theta)]])
	atomCord = np.dot(Ry,atomCord)
	return atomCord.T
	
def rotateAboutXAxis(atomCord,theta):
	atomCord = atomCord.T
	Rx = np.array([[1,0,0], [0, np.cos(theta), -np.sin(theta)], [0, np.sin(theta), np.cos(theta)]])
	atomCord = np.dot(Rx,atomCord)
	return atomCord.T
	
def findAnglebtwnVecs(n1,n2):
#finds angles between two planes given their normals
#really this function will find the angle between two vectors as well
	# for i in range(0,3):
		# if np.isnan(n1[i]):
			# n1[i] = 0
		# if np.isnan(n2[i]):
			# n2[i] = 0
	#print n1,n2
	#theta = np.absolute(np.dot(n1,n2)) / (np.linalg.norm(n1) * np.linalg.norm(n2))
	theta = np.dot(n1,n2) / (np.linalg.norm(n1) * np.linalg.norm(n2))
	theta = np.arccos(theta)
	return theta
	
def findClosestZn(x,y,z,atomCordS,atomListS):
#this function returns the two closest Zn atom coordinates given an x,y,z,
#which is presumably the coordinates of an O atom. If the closest Zn atom is far away
#because the O atom is on an edge, will return [100 100 100]
	numAtoms = len(atomListS)
	Zn1dist = 100 #obviously min distance will be much less than 100
	Zn2dist = 100
	Zn1 = 100*np.ones([1,3])
	Zn2 = 100*np.ones([1,3]) #these are just placeholders
	for i in range(0,numAtoms):
		if atomListS[i] == "Zn": #only want to calc distance if it's a Zn atom
			dx = x - atomCordS[i,0]
			dy = y - atomCordS[i,1]
			dz = z - atomCordS[i,2]
			thisDist = np.sqrt(dx**2+dy**2+dz**2)
			if thisDist < Zn1dist :
				Zn2dist = Zn1dist
				Zn1dist = thisDist
				Zn2 = Zn1
				Zn1 = np.array([atomCordS[i,0],atomCordS[i,1],atomCordS[i,2]])
			elif thisDist >= Zn1dist and thisDist < Zn2dist:
				Zn2dist = thisDist
				Zn2 = np.array([atomCordS[i,0],atomCordS[i,1],atomCordS[i,2]])
	#now check to see that distances are reasonable. if unreasonable, return placeholders
	cutOff = 4
	if Zn2dist > cutOff:
		Zn2dist = 100
		Zn2 = 100*np.ones([1,3])
	if Zn1dist > cutOff:
		Zn1dist = 100
		Zn1 = 100*np.ones([1,3]) 
	return Zn1,Zn2

def calcLinkerCenter(atomCordL,atomListL):
	numAtoms = len(atomListL)
	numInCenter = 0
	x = 0
	y = 0 
	z = 0
	for i in range(0,numAtoms):
		if atomListL[i] == "C" or atomListL[i] == "N":
			x = x + atomCordL[i,0]
			y = y + atomCordL[i,1]
			z = z + atomCordL[i,2]
			numInCenter = numInCenter + 1
	x = x/numInCenter
	y = y/numInCenter
	z = z/numInCenter
	return x,y,z

def translateCord(atomCord,tx,ty,tz):
#performs a translation of all coordinates
	T = np.array([[1, 0, 0, tx],[0,1,0,ty],[0, 0, 1, tz,],[0,0,0,1]])
	numAtoms = atomCord.shape[0]
	pad = np.ones([1,numAtoms])
	atomCord = np.append(atomCord.T,pad,axis=0)
	atomCord = np.dot(T,atomCord)
	atomCord = atomCord[0:3,:]
	return atomCord.T

	
	
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

def reduceToUnitCell(atomCord,atomList,cell_params,supercellL=0,supercellU=1):
#first convert from cartesian to fractional coordinates
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
	##go through and find which atoms to keep
	numAtoms = len(atomList)
	keep = []
	for i in range(0,numAtoms):
		xcheck = atomCordFrac[i,0] >= supercellL and atomCordFrac[i,0] < supercellU
		ycheck = atomCordFrac[i,1] >= supercellL and atomCordFrac[i,1] < supercellU
		zcheck = atomCordFrac[i,2] >= supercellL and atomCordFrac[i,2] < supercellU
		if xcheck and ycheck and zcheck:
			keep.append(i)
			
	numKeptAtoms = len(keep)
	newCord = np.empty((numKeptAtoms,3))
	#newfCord = np.empty((numKeptAtoms,3))
	newList = []
	j = 0 #index for newCord
	for i in keep:
		newCord[j,:] = atomCord[i,:]
		#newfCord[j,:] = atomCordFrac[i,:]
		j = j + 1
		newList.append(atomList[i])
	#writexyzFile(newfCord,newList,"frac.xyz")
	return  newCord, newList
		
	