#!/bin/bash

if test -z "$1"
then
  REGION=europe
else
  REGION=$1
fi

source $REGION/extents.txt

python3 ../scripts/plot_image_map.py $EXTENTS --intermediate \
	--vmin -5.0 --vmax 5.0 \
	--cm Spectral \
	-i $REGION/joint_mean_tide.txt

