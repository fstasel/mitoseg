@echo off

REM Do not specify src, dst, and settings-file arguments manually, instead:
REM 1. Create a directory (e.g. Desktop\data)
REM 2. Put input slice images into "data\src"
REM 3. Create an empty output directory "data\output"
REM 4. If needed, put custom settings into "data\settings.yaml" file
REM 5. Specify this directory in line 16 below as:
REM
REM    -v C:\path\to\your\data\directory:/app/data/
REM
REM Feel free to specify other arguments (zrange, psize, phase, etc.)
REM in line 17 below:

docker run --rm ^
    -v C:\Users\mitoseg\Desktop\data:/app/data ^
    mitoseg --zrange 30 100 --psize 2.0 dataset_slice%04d.tif
