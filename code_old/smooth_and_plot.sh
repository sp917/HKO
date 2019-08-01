#!/bin/bash


var=$1
passes=$2
mu=$3
alt=$4

printf -v mu "%g" $mu

if [ -z ${5+x} ]; then 	
	python3 main.py $var $passes $mu $alt
	cp ../plots/Smooth${passes}_${var}.png ../../../mnt/c/Users/hko/Desktop/plots/
	cp ../plots/Spectrum_smooth${passes}_${var}.png ../../../mnt/c/Users/hko/Desktop/plots/
	cp ../plots/Turbulence${alt}_added_spectra_smooth${passes}_mu${mu}_${var}.png ../../../mnt/c/Users/hko/Desktop/plots/
	cp ../plots/Turbulence${alt}_added_smooth${passes}_mu${mu}_${var}.png ../../../mnt/c/Users/hko/Desktop/plots/
else
	python3 main.py $1 $2 $3 $4
	cp ../plots/Smooth${2}_${1}_height_${3}.png ../../../mnt/c/Users/hko/Desktop/
	cp ../plots/Spectrum_smooth${2}_${1}_height${3}.png ../../../mnt/c/Users/hko/Desktop/ 
fi
echo "Done!"
