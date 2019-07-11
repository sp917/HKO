#!/bin/bash

if [ $# -eq 13 ]; then
	variable=$1
	year=$2	
	day=$3
	month=$4
	hour=$5
	layer=$6
	ow1=$7
	ow2=$8
	alt=$9
	a11=${10}
	a12=${11}
	a21=${12}
	a22=${13}
else
	variable=
	year=
	month=
	day=
	hour=
	layer=
	ow1=
	ow2=
	alt=
	a11=
	a12=
	a21=
	a22=
	echo "Plot which variable? Options are U, V, ke, speed, U10, V10, speeed10, ke10."
	read variable
	echo "For which year?"
	read year
	echo "Which day?"
	read day
	echo "Which month?"
	read month
	echo "Which hour?"
	read hour
	echo "Which layer?"
	read layer
	echo "Overwrite interpolated field data? Type 1 for yes."
	read ow1
	echo "Overwrite spectrum data? Type 1 for yes."
	read ow2
	echo "Use ncl routine for calculating spectra? Type 1 for yes."
	read alt
	echo "Enter a11:"
	read a11
	echo "Enter a12:"
	read a12
	echo "Enter a21:"
	read a21
	echo "Enter a22:"
	read a22

fi

if [ $ow1 -eq 1 ]; then
	conda activate ncl_stable	
	for t in 00 01 12
	do
		echo "Interpolating data onto vertical levels..."
		ncl yy=$year mm=$month dd=$day hh=$hour tt=$t ./vert_interp.ncl
		echo "Done."
	done
fi

conda activate base
if [ $ow2 -eq 1 ]; then
	for t in 00 01 12
	do
		echo "Beginning python routine to calculate spectrum..."
		python3 readdata.py $variable $year $month $day $hour $t
		echo "Done."
	done
fi


python3 plotall.py $variable $year $month $day $hour $layer $a11 $a12 $a21 $a22 

printf -v year "%04d" $year
printf -v month "%02d" $month
printf -v day "%02d" $day
printf -v hour "%02d" $hour

cp ./plots/spectra_${variable}_layer_${layer}_${year}-${day}-${month}_${hour}:00:00.png ../../mnt/c/Users/hko/Desktop/

cp ./plots/alt_spectra_${variable}_layer_${layer}_${year}-${day}-${month}_${hour}:00:00.png ../../mnt/c/Users/hko/Desktop/

echo "Done!"

