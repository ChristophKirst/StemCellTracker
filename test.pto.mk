
# makefile for panorama stitching, created by hugin using the new makefilelib

# Tool configuration
NONA=nona
PTSTITCHER=PTStitcher
PTMENDER=PTmender
PTBLENDER=PTblender
PTMASKER=PTmasker
PTROLLER=PTroller
ENBLEND=enblend
ENFUSE=enfuse
SMARTBLEND=smartblend.exe
HDRMERGE=hugin_hdrmerge
RM=rm
EXIFTOOL=exiftool

# Project parameters
HUGIN_PROJECTION=2
HUGIN_HFOV=3
HUGIN_WIDTH=1536
HUGIN_HEIGHT=1024

# options for the programs
NONA_LDR_REMAPPED_COMP=-z LZW
NONA_OPTS=
ENBLEND_OPTS= -f943x768+512+256
ENBLEND_LDR_COMP=--compression=LZW
ENBLEND_EXPOSURE_COMP=--compression=LZW
ENBLEND_HDR_COMP=
HDRMERGE_OPTS=
ENFUSE_OPTS=
EXIFTOOL_COPY_ARGS=-ImageDescription -Make -Model -Artist -WhitePoint -Copyright -GPS:all -DateTimeOriginal -CreateDate -UserComment -ColorSpace -OwnerName -SerialNumber
EXIFTOOL_INFO_ARGS='-Software=Hugin 2013.0.0.4692917e7a55' '-UserComment<$${UserComment}&\#xa;Projection: Equirectangular (2)&\#xa;FOV: 3 x 2&\#xa;Ev: 0.00' -f

# the output panorama
LDR_REMAPPED_PREFIX=test
LDR_REMAPPED_PREFIX_SHELL=test
HDR_STACK_REMAPPED_PREFIX=test_hdr_
HDR_STACK_REMAPPED_PREFIX_SHELL=test_hdr_
LDR_EXPOSURE_REMAPPED_PREFIX=test_exposure_layers_
LDR_EXPOSURE_REMAPPED_PREFIX_SHELL=test_exposure_layers_
PROJECT_FILE=/home/ckirst/Science/Simulation/Matlab/StemCell3D/test.pto
PROJECT_FILE_SHELL=/home/ckirst/Science/Simulation/Matlab/StemCell3D/test.pto
LDR_BLENDED=test.tif
LDR_BLENDED_SHELL=test.tif
LDR_STACKED_BLENDED=test_fused.tif
LDR_STACKED_BLENDED_SHELL=test_fused.tif
LDR_EXPOSURE_LAYERS_FUSED=test_blended_fused.tif
LDR_EXPOSURE_LAYERS_FUSED_SHELL=test_blended_fused.tif
HDR_BLENDED=test_hdr.exr
HDR_BLENDED_SHELL=test_hdr.exr

# first input image
INPUT_IMAGE_1=/home/ckirst/Science/Simulation/Matlab/StemCell3D/test0001.tif
INPUT_IMAGE_1_SHELL=/home/ckirst/Science/Simulation/Matlab/StemCell3D/test0001.tif

# all input images
INPUT_IMAGES=/home/ckirst/Science/Simulation/Matlab/StemCell3D/test0001.tif\
/home/ckirst/Science/Simulation/Matlab/StemCell3D/test0002.tif\
/home/ckirst/Science/Simulation/Matlab/StemCell3D/test0003.tif\
/home/ckirst/Science/Simulation/Matlab/StemCell3D/test0004.tif
INPUT_IMAGES_SHELL=/home/ckirst/Science/Simulation/Matlab/StemCell3D/test0001.tif\
/home/ckirst/Science/Simulation/Matlab/StemCell3D/test0002.tif\
/home/ckirst/Science/Simulation/Matlab/StemCell3D/test0003.tif\
/home/ckirst/Science/Simulation/Matlab/StemCell3D/test0004.tif

# remapped images
LDR_LAYERS=test0000.tif\
test0001.tif\
test0002.tif\
test0003.tif
LDR_LAYERS_SHELL=test0000.tif\
test0001.tif\
test0002.tif\
test0003.tif

# remapped images (hdr)
HDR_LAYERS=test_hdr_0000.exr\
test_hdr_0001.exr\
test_hdr_0002.exr\
test_hdr_0003.exr
HDR_LAYERS_SHELL=test_hdr_0000.exr\
test_hdr_0001.exr\
test_hdr_0002.exr\
test_hdr_0003.exr

# remapped maxval images
HDR_LAYERS_WEIGHTS=test_hdr_0000_gray.pgm\
test_hdr_0001_gray.pgm\
test_hdr_0002_gray.pgm\
test_hdr_0003_gray.pgm
HDR_LAYERS_WEIGHTS_SHELL=test_hdr_0000_gray.pgm\
test_hdr_0001_gray.pgm\
test_hdr_0002_gray.pgm\
test_hdr_0003_gray.pgm

# stacked hdr images
HDR_STACK_0=test_stack_hdr_0000.exr
HDR_STACK_0_SHELL=test_stack_hdr_0000.exr
HDR_STACK_0_INPUT=test_hdr_0000.exr
HDR_STACK_0_INPUT_SHELL=test_hdr_0000.exr
HDR_STACK_1=test_stack_hdr_0001.exr
HDR_STACK_1_SHELL=test_stack_hdr_0001.exr
HDR_STACK_1_INPUT=test_hdr_0001.exr
HDR_STACK_1_INPUT_SHELL=test_hdr_0001.exr
HDR_STACK_2=test_stack_hdr_0002.exr
HDR_STACK_2_SHELL=test_stack_hdr_0002.exr
HDR_STACK_2_INPUT=test_hdr_0002.exr
HDR_STACK_2_INPUT_SHELL=test_hdr_0002.exr
HDR_STACK_3=test_stack_hdr_0003.exr
HDR_STACK_3_SHELL=test_stack_hdr_0003.exr
HDR_STACK_3_INPUT=test_hdr_0003.exr
HDR_STACK_3_INPUT_SHELL=test_hdr_0003.exr
HDR_STACKS_NUMBERS=0 1 2 3 
HDR_STACKS=$(HDR_STACK_0) $(HDR_STACK_1) $(HDR_STACK_2) $(HDR_STACK_3) 
HDR_STACKS_SHELL=$(HDR_STACK_0_SHELL) $(HDR_STACK_1_SHELL) $(HDR_STACK_2_SHELL) $(HDR_STACK_3_SHELL) 

# number of image sets with similar exposure
LDR_EXPOSURE_LAYER_0=test_exposure_0000.tif
LDR_EXPOSURE_LAYER_0_SHELL=test_exposure_0000.tif
LDR_EXPOSURE_LAYER_0_INPUT=test_exposure_layers_0000.tif\
test_exposure_layers_0001.tif\
test_exposure_layers_0002.tif\
test_exposure_layers_0003.tif
LDR_EXPOSURE_LAYER_0_INPUT_SHELL=test_exposure_layers_0000.tif\
test_exposure_layers_0001.tif\
test_exposure_layers_0002.tif\
test_exposure_layers_0003.tif
LDR_EXPOSURE_LAYER_0_INPUT_PTMENDER=test0000.tif\
test0001.tif\
test0002.tif\
test0003.tif
LDR_EXPOSURE_LAYER_0_INPUT_PTMENDER_SHELL=test0000.tif\
test0001.tif\
test0002.tif\
test0003.tif
LDR_EXPOSURE_LAYER_0_EXPOSURE=0
LDR_EXPOSURE_LAYERS_NUMBERS=0 
LDR_EXPOSURE_LAYERS=$(LDR_EXPOSURE_LAYER_0) 
LDR_EXPOSURE_LAYERS_SHELL=$(LDR_EXPOSURE_LAYER_0_SHELL) 
LDR_EXPOSURE_LAYERS_REMAPPED=test_exposure_layers_0000.tif\
test_exposure_layers_0001.tif\
test_exposure_layers_0002.tif\
test_exposure_layers_0003.tif
LDR_EXPOSURE_LAYERS_REMAPPED_SHELL=test_exposure_layers_0000.tif\
test_exposure_layers_0001.tif\
test_exposure_layers_0002.tif\
test_exposure_layers_0003.tif

# stacked ldr images
LDR_STACK_0=test_stack_ldr_0000.tif
LDR_STACK_0_SHELL=test_stack_ldr_0000.tif
LDR_STACK_0_INPUT=test_exposure_layers_0000.tif
LDR_STACK_0_INPUT_SHELL=test_exposure_layers_0000.tif
LDR_STACK_1=test_stack_ldr_0001.tif
LDR_STACK_1_SHELL=test_stack_ldr_0001.tif
LDR_STACK_1_INPUT=test_exposure_layers_0001.tif
LDR_STACK_1_INPUT_SHELL=test_exposure_layers_0001.tif
LDR_STACK_2=test_stack_ldr_0002.tif
LDR_STACK_2_SHELL=test_stack_ldr_0002.tif
LDR_STACK_2_INPUT=test_exposure_layers_0002.tif
LDR_STACK_2_INPUT_SHELL=test_exposure_layers_0002.tif
LDR_STACK_3=test_stack_ldr_0003.tif
LDR_STACK_3_SHELL=test_stack_ldr_0003.tif
LDR_STACK_3_INPUT=test_exposure_layers_0003.tif
LDR_STACK_3_INPUT_SHELL=test_exposure_layers_0003.tif
LDR_STACKS_NUMBERS=0 1 2 3 
LDR_STACKS=$(LDR_STACK_0) $(LDR_STACK_1) $(LDR_STACK_2) $(LDR_STACK_3) 
LDR_STACKS_SHELL=$(LDR_STACK_0_SHELL) $(LDR_STACK_1_SHELL) $(LDR_STACK_2_SHELL) $(LDR_STACK_3_SHELL) 
DO_LDR_BLENDED=1
DO_LDR_EXPOSURE_LAYERS_FUSED=1

all : startStitching $(LDR_BLENDED) $(LDR_LAYERS) $(LDR_EXPOSURE_LAYERS_FUSED) 

startStitching : 
	@echo '==========================================================================='
	@echo 'Stitching panorama'
	@echo '==========================================================================='

clean : 
	@echo '==========================================================================='
	@echo 'Remove temporary files'
	@echo '==========================================================================='
	-$(RM) $(LDR_EXPOSURE_LAYERS_SHELL) $(LDR_EXPOSURE_LAYERS_REMAPPED_SHELL) 

test : 
	@echo '==========================================================================='
	@echo 'Testing programs'
	@echo '==========================================================================='
	@echo -n 'Checking nona...'
	@-$(NONA) --help > /dev/null 2>&1 && echo '[OK]' || echo '[FAILED]'
	@echo -n 'Checking enblend...'
	@-$(ENBLEND) -h > /dev/null 2>&1 && echo '[OK]' || echo '[FAILED]'
	@echo -n 'Checking enfuse...'
	@-$(ENFUSE) -h > /dev/null 2>&1 && echo '[OK]' || echo '[FAILED]'
	@echo -n 'Checking hugin_hdrmerge...'
	@-$(HDRMERGE) -h > /dev/null 2>&1 && echo '[OK]' || echo '[FAILED]'
	@echo -n 'Checking exiftool...'
	@-$(EXIFTOOL) -ver > /dev/null 2>&1 && echo '[OK]' || echo '[FAILED]'

info : 
	@echo '==========================================================================='
	@echo '***************  Panorama makefile generated by Hugin       ***************'
	@echo '==========================================================================='
	@echo 'System information'
	@echo '==========================================================================='
	@echo -n 'Operating system: '
	@-uname -o
	@echo -n 'Release: '
	@-uname -r
	@echo -n 'Kernel version: '
	@-uname -v
	@echo -n 'Machine: '
	@-uname -m
	@echo 'Disc usage'
	@-df -h
	@echo 'Memory usage'
	@-free -m
	@echo '==========================================================================='
	@echo 'Output options'
	@echo '==========================================================================='
	@echo 'Hugin Version: 2013.0.0.4692917e7a55'
	@echo 'Project file: /home/ckirst/Science/Simulation/Matlab/StemCell3D/test.pto'
	@echo 'Output prefix: test'
	@echo 'Projection: Equirectangular (2)'
	@echo 'Field of view: 3 x 2'
	@echo 'Canvas dimensions: 1536 x 1024'
	@echo 'Crop area: (512,256) - (1455,1024)'
	@echo 'Output exposure value: 0.00'
	@echo 'Output stacks minimum overlap: 0.700'
	@echo 'Output layers maximum Ev difference: 0.50'
	@echo 'Selected outputs'
	@echo 'Normal panorama'
	@echo '* Blended panorama'
	@echo '* Remapped images'
	@echo 'Exposure fusion'
	@echo '* Blended and fused panorama'
	@echo '==========================================================================='
	@echo 'Input images'
	@echo '==========================================================================='
	@echo 'Number of images in project file: 4'
	@echo 'Number of active images: 4'
	@echo 'Image 0: /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0001.tif'
	@echo 'Image 0: Size 512x512, Exposure: 0.00'
	@echo 'Image 1: /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0002.tif'
	@echo 'Image 1: Size 512x512, Exposure: 0.00'
	@echo 'Image 2: /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0003.tif'
	@echo 'Image 2: Size 512x512, Exposure: 0.00'
	@echo 'Image 3: /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0004.tif'
	@echo 'Image 3: Size 512x512, Exposure: 0.00'

# Rules for ordinary TIFF_m and hdr output

test0000.tif : /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0001.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) $(NONA_LDR_REMAPPED_COMP) -r ldr -m TIFF_m -o $(LDR_REMAPPED_PREFIX_SHELL) -i 0 $(PROJECT_FILE_SHELL)

test_hdr_0000.exr : /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0001.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) -r hdr -m EXR_m -o $(HDR_STACK_REMAPPED_PREFIX_SHELL) -i 0 $(PROJECT_FILE_SHELL)

test0001.tif : /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0002.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) $(NONA_LDR_REMAPPED_COMP) -r ldr -m TIFF_m -o $(LDR_REMAPPED_PREFIX_SHELL) -i 1 $(PROJECT_FILE_SHELL)

test_hdr_0001.exr : /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0002.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) -r hdr -m EXR_m -o $(HDR_STACK_REMAPPED_PREFIX_SHELL) -i 1 $(PROJECT_FILE_SHELL)

test0002.tif : /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0003.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) $(NONA_LDR_REMAPPED_COMP) -r ldr -m TIFF_m -o $(LDR_REMAPPED_PREFIX_SHELL) -i 2 $(PROJECT_FILE_SHELL)

test_hdr_0002.exr : /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0003.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) -r hdr -m EXR_m -o $(HDR_STACK_REMAPPED_PREFIX_SHELL) -i 2 $(PROJECT_FILE_SHELL)

test0003.tif : /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0004.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) $(NONA_LDR_REMAPPED_COMP) -r ldr -m TIFF_m -o $(LDR_REMAPPED_PREFIX_SHELL) -i 3 $(PROJECT_FILE_SHELL)

test_hdr_0003.exr : /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0004.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) -r hdr -m EXR_m -o $(HDR_STACK_REMAPPED_PREFIX_SHELL) -i 3 $(PROJECT_FILE_SHELL)

# Rules for exposure layer output

test_exposure_layers_0000.tif : /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0001.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) $(NONA_LDR_REMAPPED_COMP) -r ldr -e 0 -m TIFF_m -o $(LDR_EXPOSURE_REMAPPED_PREFIX_SHELL) -i 0 $(PROJECT_FILE_SHELL)

test_exposure_layers_0001.tif : /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0002.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) $(NONA_LDR_REMAPPED_COMP) -r ldr -e 0 -m TIFF_m -o $(LDR_EXPOSURE_REMAPPED_PREFIX_SHELL) -i 1 $(PROJECT_FILE_SHELL)

test_exposure_layers_0002.tif : /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0003.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) $(NONA_LDR_REMAPPED_COMP) -r ldr -e 0 -m TIFF_m -o $(LDR_EXPOSURE_REMAPPED_PREFIX_SHELL) -i 2 $(PROJECT_FILE_SHELL)

test_exposure_layers_0003.tif : /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0004.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) $(NONA_LDR_REMAPPED_COMP) -r ldr -e 0 -m TIFF_m -o $(LDR_EXPOSURE_REMAPPED_PREFIX_SHELL) -i 3 $(PROJECT_FILE_SHELL)

# Rules for LDR and HDR stack merging, a rule for each stack

$(LDR_STACK_0) : $(LDR_STACK_0_INPUT) 
	$(ENFUSE) $(ENFUSE_OPTS) -o $(LDR_STACK_0_SHELL) -- $(LDR_STACK_0_INPUT_SHELL)
	-$(EXIFTOOL) -overwrite_original_in_place -TagsFromFile $(INPUT_IMAGE_1_SHELL) $(EXIFTOOL_COPY_ARGS) $(LDR_STACK_0_SHELL)

$(HDR_STACK_0) : $(HDR_STACK_0_INPUT) 
	$(HDRMERGE) $(HDRMERGE_OPTS) -o $(HDR_STACK_0_SHELL) -- $(HDR_STACK_0_INPUT_SHELL)

$(LDR_STACK_1) : $(LDR_STACK_1_INPUT) 
	$(ENFUSE) $(ENFUSE_OPTS) -o $(LDR_STACK_1_SHELL) -- $(LDR_STACK_1_INPUT_SHELL)
	-$(EXIFTOOL) -overwrite_original_in_place -TagsFromFile $(INPUT_IMAGE_1_SHELL) $(EXIFTOOL_COPY_ARGS) $(LDR_STACK_1_SHELL)

$(HDR_STACK_1) : $(HDR_STACK_1_INPUT) 
	$(HDRMERGE) $(HDRMERGE_OPTS) -o $(HDR_STACK_1_SHELL) -- $(HDR_STACK_1_INPUT_SHELL)

$(LDR_STACK_2) : $(LDR_STACK_2_INPUT) 
	$(ENFUSE) $(ENFUSE_OPTS) -o $(LDR_STACK_2_SHELL) -- $(LDR_STACK_2_INPUT_SHELL)
	-$(EXIFTOOL) -overwrite_original_in_place -TagsFromFile $(INPUT_IMAGE_1_SHELL) $(EXIFTOOL_COPY_ARGS) $(LDR_STACK_2_SHELL)

$(HDR_STACK_2) : $(HDR_STACK_2_INPUT) 
	$(HDRMERGE) $(HDRMERGE_OPTS) -o $(HDR_STACK_2_SHELL) -- $(HDR_STACK_2_INPUT_SHELL)

$(LDR_STACK_3) : $(LDR_STACK_3_INPUT) 
	$(ENFUSE) $(ENFUSE_OPTS) -o $(LDR_STACK_3_SHELL) -- $(LDR_STACK_3_INPUT_SHELL)
	-$(EXIFTOOL) -overwrite_original_in_place -TagsFromFile $(INPUT_IMAGE_1_SHELL) $(EXIFTOOL_COPY_ARGS) $(LDR_STACK_3_SHELL)

$(HDR_STACK_3) : $(HDR_STACK_3_INPUT) 
	$(HDRMERGE) $(HDRMERGE_OPTS) -o $(HDR_STACK_3_SHELL) -- $(HDR_STACK_3_INPUT_SHELL)

$(LDR_BLENDED) : $(LDR_LAYERS) 
	$(ENBLEND) $(ENBLEND_LDR_COMP) $(ENBLEND_OPTS) -o $(LDR_BLENDED_SHELL) -- $(LDR_LAYERS_SHELL)
	-$(EXIFTOOL) -E -overwrite_original_in_place -TagsFromFile $(INPUT_IMAGE_1_SHELL) $(EXIFTOOL_COPY_ARGS) $(EXIFTOOL_INFO_ARGS) $(LDR_BLENDED_SHELL)

$(LDR_EXPOSURE_LAYER_0) : $(LDR_EXPOSURE_LAYER_0_INPUT) 
	$(ENBLEND) $(ENBLEND_EXPOSURE_COMP) $(ENBLEND_OPTS) -o $(LDR_EXPOSURE_LAYER_0_SHELL) -- $(LDR_EXPOSURE_LAYER_0_INPUT_SHELL)
	-$(EXIFTOOL) -overwrite_original_in_place -TagsFromFile $(INPUT_IMAGE_1_SHELL) $(EXIFTOOL_COPY_ARGS) $(LDR_EXPOSURE_LAYER_0_SHELL)

$(LDR_STACKED_BLENDED) : $(LDR_STACKS) 
	$(ENBLEND) $(ENBLEND_LDR_COMP) $(ENBLEND_OPTS) -o $(LDR_STACKED_BLENDED_SHELL) -- $(LDR_STACKS_SHELL)
	-$(EXIFTOOL) -E -overwrite_original_in_place -TagsFromFile $(INPUT_IMAGE_1_SHELL) $(EXIFTOOL_COPY_ARGS) $(EXIFTOOL_INFO_ARGS) $(LDR_STACKED_BLENDED_SHELL)

$(LDR_EXPOSURE_LAYERS_FUSED) : $(LDR_EXPOSURE_LAYERS) 
	$(ENFUSE) $(ENBLEND_LDR_COMP) $(ENFUSE_OPTS) -o $(LDR_EXPOSURE_LAYERS_FUSED_SHELL) -- $(LDR_EXPOSURE_LAYERS_SHELL)
	-$(EXIFTOOL) -E -overwrite_original_in_place -TagsFromFile $(INPUT_IMAGE_1_SHELL) $(EXIFTOOL_COPY_ARGS) $(EXIFTOOL_INFO_ARGS) $(LDR_EXPOSURE_LAYERS_FUSED_SHELL)

$(HDR_BLENDED) : $(HDR_STACKS) 
	$(ENBLEND) $(ENBLEND_HDR_COMP) $(ENBLEND_OPTS) -o $(HDR_BLENDED_SHELL) -- $(HDR_STACKS_SHELL)

$(LDR_REMAPPED_PREFIX)_multilayer.tif : $(LDR_LAYERS) 
	tiffcp $(LDR_LAYERS_SHELL) $(LDR_REMAPPED_PREFIX_SHELL)_multilayer.tif

$(LDR_REMAPPED_PREFIX)_fused_multilayer.tif : $(LDR_STACKS) $(LDR_EXPOSURE_LAYERS) 
	tiffcp $(LDR_STACKS_SHELL) $(LDR_EXPOSURE_LAYERS_SHELL) $(LDR_REMAPPED_PREFIX_SHELL)_fused_multilayer.tif

$(LDR_REMAPPED_PREFIX)_multilayer.psd : $(LDR_LAYERS) 
	PTtiff2psd -o $(LDR_REMAPPED_PREFIX_SHELL)_multilayer.psd $(LDR_LAYERS_SHELL)

$(LDR_REMAPPED_PREFIX)_fused_multilayer.psd : $(LDR_STACKS) $(LDR_EXPOSURE_LAYERS) 
	PTtiff2psd -o $(LDR_REMAPPED_PREFIX_SHELL)_fused_multilayer.psd $(LDR_STACKS_SHELL)$(LDR_EXPOSURE_LAYERS_SHELL)
