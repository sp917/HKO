#!/bin/bash

if [ $# -eq 8 ]; then
	variable=$1
	year=$2	
	day=$3
	month=$4
	hour=$5
	layer=$6
	a1=$7
	a2=$8
else
	variable=
	year=
	month=
	day=
	hour=
	layer=
	a=
	echo "Plot which variable? Options are U, V, ke, speed."
	read variable
	echo "For which year?"
	read year
	echo "Which day?"
	read day
	echo "Which month?"
	read month
	echo "Which hour?"
	read hour
	echo "On which layer?"
	read layer
	echo "Enter a value for a"
	read a
fi

printf -v year "%04d" $year
printf -v month "%02d" $month
printf -v day "%02d" $day
printf -v hour "%02d" $hour

echo "Plotting spectra..."

python3 plotall.py $variable $year $month $day $hour $layer $a1 $a2

cp ./plots/spectra_${variable}_layer_${layer}_${year}-${day}-${month}_${hour}:00:00.png ../../mnt/c/Users/hko/Desktop/

echo "Done!"
