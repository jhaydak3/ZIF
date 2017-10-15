#!/usr/bin/env python 
import numpy as np
import sys



def main(file_to_read,mol_track_ndx,SAMPLE_FRAME_RATE):
	fh=open(file_to_read,"r")
	fh_new=open(file_to_read[:-4] + "_traj.xyz","w")
	#fh_yz=open(file_to_read[:-4] + "_yz.txt","w")
	#fh_xz=open(file_to_read[:-4] + "_xz.txt","w")
	fh_xyz=open(file_to_read[:-4] + "_xyz.txt","w")
	#read the first three lines, they are junk
	_ = fh.readline()
	_ = fh.readline()
	_ = fh.readline()
	breakWhileLoop = False
	j = 1
	while not breakWhileLoop:
		if j%SAMPLE_FRAME_RATE == 0:
			print j
			thisLine = fh.readline()
			if not thisLine:
				break
			_, numCOM = thisLine.split()
			numCOM = int(numCOM)
			fh_new.write(str(numCOM) + "\n")
			fh_new.write("For Whom The Bell Tolls \n") ## doesnt matter
			for i in range(0,numCOM):
				thisLine = fh.readline()
	   
				if not thisLine: ##end of file
					breakWhileLoop = True
					break 
		
				fh_new.write(thisLine)
				if i == mol_track_ndx-1:
					_, x, y, z = thisLine.split()
					#fh_yz.write(y + " " + z + "\n")
					#fh_xz.write(x + " " + z + "\n")
					fh_xyz.write(x + " " + y + " "+ z+"\n")
		else:
			thisLine = fh.readline()
			if not thisLine:
				break
			_, numCOM = thisLine.split()
			numCOM = int(numCOM)
			for k in range(0,numCOM):
				_ = fh.readline()
				if not thisLine: ##end of file
					breakWhileLoop = True
					break 
			
		j = j+1
	fh.close()
	fh_new.close()
	#fh_yz.close()
	#fh_xz.close()
	fh_xyz.close()    
     
 
if __name__ == "__main__":
	#INPUTS
	mol_track_ndx = 4 #which molecule to track?
	file_to_read = "/nv/pj1/jhaydak3/FORCE_FIELD_MAKER/lmp_files/BzIm-ZIF-7-CO2-run1/MOLECULE_msd.OUT"
	SAMPLE_FRAME_RATE=50 #only keep every 10 frames
	main(file_to_read,mol_track_ndx,SAMPLE_FRAME_RATE)
