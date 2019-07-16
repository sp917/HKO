#!/bin/bash

if [ $# -eq 10 ]; then
	variable=$1
	year=$2	
	day=$3
	month=$4
	hour=$5
	layer=$6
	ny=$7
	nx=$8
	ow=$9
	filt=${10}
else
	variable=
	year=
	month=
	day=
	hour=
	layer=
	ny=
	nx=
	ow=
	filt=
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
	echo "Enter ny:"
	read ny
	echo "Enter nx:"
	read nx
	echo "Overwrite filter data? Type 1 for yes."
	read ow
	echo "Which filter to use? Type 1 or 2."
	read filt
fi

printf -v year "%04d" $year
printf -v month "%02d" $month
printf -v day "%02d" $day
printf -v hour "%02d" $hour

echo $ow

if [ $ow -eq 1 ]; then
	echo "Filtering"
	python3 ./filter.py $variable $year $month $day $hour $ny $nx $filt
fi

echo "Plotting fields..."

python3 ./filterplot.py $variable $year $month $day $hour $layer $ny $nx $filt

cp ./plots/Filtered${filt}_${ny}_${nx}_${variable}_layer_${layer}_${year}-${day}-${month}_${hour}:00:00.png ../../mnt/c/Users/hko/Desktop/
cp ./plots/Filtered${filt}_${ny}_${nx}_${variable}_${year}-${day}-${month}_${hour}:00:00.png ../../mnt/c/Users/hko/Desktop/


echo "Done!"
