#!/bin/bash



if [ -z ${3+x} ]; then 	
	python3 main.py $1 $2 
	cp ../plots/Smooth${2}_${1}.png ../../../mnt/c/Users/hko/Desktop/
	cp ../plots/Spectrum_smooth${2}_$1.png ../../../mnt/c/Users/hko/Desktop/
else
	python3 main.py $1 $2 $3
	cp ../plots/Smooth${2}_${1}_height_${3}.png ../../../mnt/c/Users/hko/Desktop/
	cp ../plots/Spectrum_smooth${2}_${1}_height${3}.png ../../../mnt/c/Users/hko/Desktop/ 
fi
echo "Done!"
