#!/bin/bash

# Do not specify src, dst, and settingsFile arguments manually, instead:
# 1. Create a directory (e.g. Desktop/data)
# 2. Put input slice images into "data/src"
# 3. Create an empty output directory "data/output"
# 4. If needed, put custom settings into "data/settings.yaml" file
# 5. Specify this directory in line 16 below as:
#
#    -v /path/to/your/data/directory:/app/data/ \
#
# Feel free to specify other arguments (psize, zrange, phase, etc.)
# in line 17 below:

docker run --rm \
	-v /home/mitoseg/Desktop/data:/app/data/ \
	mitoseg -psize 2.0 -zrange 30 100 dataset_slice%04d.tif
