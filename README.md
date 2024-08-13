# MitoSeg - Mitochondria Segmentation Tool

Based on the paper:

> Tasel, S.F., Mumcuoglu, E.U., Hassanpour, R.Z. and Perkins, G., 2016. A validated active contour method driven by parabolic arc model for detection and segmentation of mitochondria. Journal of structural biology, 194(3), pp.253-271.

Please cite this paper if you use MitoSeg in your research.

## Setup

MitoSeg uses the following libraries:

- [OpenCV](https://opencv.org/)
- [Boost](https://www.boost.org/)
- [yaml-cpp](https://github.com/jbeder/yaml-cpp)

To compile and use MitoSeg, development files for OpenCV4, Boost, and yaml-cpp must be installed. For GNU/Linux distributions based on Debian (e.g., Ubuntu):

    apt install libopencv-dev libboost-dev libboost-program-options-dev libyaml-cpp-dev

Please refer to respective user guides for other GNU/Linux distributions or operating systems.

MitoSeg can be compiled with:

    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make -j $(nproc)

## Usage

### Synopsis

    mitoseg [OPTION] --zrange # # --psize # FILENAME_PATTERN

### Mandatory Arguments

| Argument | Description |
| -- | -- |
| `--zrange <start slice #> <end slice #>` | Slice z-range start / end. |
| `--psize <pixel size>` | Pixel size in nm/px. |
| `FILENAME_PATTERN` | Filename of slice images and slice number tag. `%d` indicates the position of the slice number tag. For example, `mito%d.bmp` generates file names such as `mito0.bmp`, `mito1.bmp`, etc. and `mito%03d.bmp` generates file names such as `mito000.bmp`, `mito001.bmp`, etc. |


### Optional Arguments

| Argument | Description |
| -- | -- |
| `--roi <left> <top> <width> <height>` | Specifies the left and top pixel coordinates and the width and height of the bounding rectangle of the interested region of segmentation. If not used, ROI will be determined automatically. |
| `--src <directory>` | Specifies the path of the dataset image files. |
| `--dst <directory>` | Specifies the directory where the outputs (intermediate data and image files, final segmentation image files, an [IMOD](https://bio3d.colorado.edu/imod/) model file (.mod), and a 3D mesh file (.ply)) are stored. |
| `--valid <validity threshold>` | Specifies the validity threshold that is explained in the paper. Must be between 0 and 1 (default: 0.75). |
| `--thick <z-thickness>` | Specifies the snake z-thickness. The given value must be between 5 and 500 (default: 20). Can be set to `full` to use the whole z-range specified by the `--zrange` argument. Snake thickness affects the segmentation accuracy; thicker snakes provide continuous smooth segmentation along the z-range, but can also produce more false negatives. |
| `--phase <phase #>` | Specifies a specific phase to be executed (1, 2, or 3). The `--valid` and `--thick` arguments affect the 2nd and 3rd phases, respectively. Hence, if the same portion of the dataset is to be segmented by changing only these arguments, the `--phase` argument can be used to restart a particular phase to avoid redundant computation manually. If not used, MitoSeg will execute all phases in order. |
| `--cores` | Specifies the number of CPU cores that work in collaboration to speed up the process. |
| `--settingsFile <file path>` | Specifies the path of a YAML file containing custom segmentation variables. |

## Examples

- Process the files dataset_slice0030.tif to dataset_slice0100.tif assuming that pixel size is 2.0nm:

      ./mitoseg --zrange 30 100 --psize 2.0 dataset_slice%04d.tif

- Process the files mito40.bmp to mito120.bmp from the slices directory, assuming that pixel size is 1.1nm:
   
      ./mitoseg --src slices --zrange 40 120 --psize 1.1 mito%d.bmp

- Load custom settings from settings.yaml, process the files slice0020.bmp to slice0080.bmp, assuming the pixel size is 2.1nm, and save the results in the outputs directory:

      ./mitoseg --settingsFile settings.yaml --dst outputs --zrange 20 80 --psize 2.1 slice00%d.bmp

## Running MitoSeg as a Docker Application
MitoSeg can also be used within a [Docker](https://www.docker.com/) environment, which eliminates the need to prepare the building environment on your computer, and allows MitoSeg to run on various operating systems easily. A sample Dockerfile is provided under the `docker` directory, which can be used for building a Docker image by executing the following command within the `docker` directory:

    docker build -t mitoseg .

After the image is built, it can be run by executing the provided `docker-mitoseg.sh` (for Linux/Mac) or `docker-mitoseg.cmd` (for Windows) scripts. These scripts include instructions and guidelines, and can be modified according to specific needs.
