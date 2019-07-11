#!/bin/bash

if [ $# -eq 10 ]; then
	variable=$1
	year=$2	
	day=$3
	month=$4
	hour=$5
	t=$6
	layer=$7
	overwrite1=$8
	overwrite2=$9
	a=${10}
else
	variable=
	year=
	month=
	day=
	hour=
	t=
	layer=
	overwrite1=
	overwrite2=
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
	echo "What time?"
	read t
	echo "On which layer?"
	read layer
	echo "Overwrite interpolated UV data? Type 1 for yes."
	read overwrite1
	echo "Overwrite spectrum data? Type 1 for yes."
	read overwrite2
	echo "Enter a value for a"
	read a
fi

printf -v year "%04d" $year
printf -v month "%02d" $month
printf -v day "%02d" $day
printf -v hour "%02d" $hour
printf -v t "%02d" $t

#varno=
#if [ $variable == "U" ]; then
#	varno=1
#elif [ $variable == "V" ]; then
#	varno=2
#elif [ $variable == "ke" ]; then
#	varno=3
#elif [ $variable == "speed" ]; then
#	varno=4
#else
#        variable=U	
#	varno=1
#fi


if [ $overwrite1 -eq 1 ]; then
	conda activate ncl_stable	
	echo "Interpolating data onto vertical levels..."
	ncl yy=$year mm=$month dd=$day hh=$hour tt=$t ./vert_interp.ncl
	echo "Done."
fi

conda activate base
if [ $overwrite2 -eq 1 ]; then
	echo "Beginning python routine to calculate spectrum..."
	python3 readdata.py $variable $year $month $day $hour $t $layer 
	echo "Done."
fi

echo "Plotting spectrum..."

# ncl var=$varno yy=$year mm=$month dd=$day hh=$hour tt=$t ll=$layer aa=$a ./plot_spectrum.ncl

python3 plot_spectrum.py $variable $year $month $day $hour $t $layer $a

cp ./plots/spectrum_${variable}_layer_${layer}_${year}-${day}-${month}_${hour}:00:00.t${t}.png ../../mnt/c/Users/hko/Desktop/

echo "Done!"
