#!/usr/bin/env bash
S_FOLDER=POC_25
FUN_NAME="JonTest1.py"
LINKER_NAME="part_Im.xyz"
HELP_FILE_NAME="readFiles.py"

XYZFILES=$(ls $S_FOLDER/*xyz)


counter=0
for THISFILE in $XYZFILES
do
	echo $THISFILE
	python JonTest1.py $THISFILE
	((counter++))
	echo Completed making unit cells for $counter files.
	
	if [ $counter -gt 25 ]
	then 
			break
	fi
done

UNIT_CELL_FILES=$(ls $S_FOLDER/*unitcell.xyz)
mv $UNIT_CELL_FILES $PWD
mkdir Zifunitcells
mv *unitcell.xyz Zifunitcells



