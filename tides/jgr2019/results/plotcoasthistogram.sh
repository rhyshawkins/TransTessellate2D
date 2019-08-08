#!/bin/bash

if test -z "$1"
then
  REGION=europe
else
  REGION=$1
fi

source $REGION/extents.txt

python3 ../scripts/plot_coast_histogram.py -H $REGION/joint_histogram_tide.txt \
	--log \
	$EXTENTS \
	--colorbar --vmax -2.0 \
	--width 8 --height 3 \
	--cmap jet \
	-c ../points/${REGION}_points.txt \
	-C ../points/${REGION}_places.txt 
