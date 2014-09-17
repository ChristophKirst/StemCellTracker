
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
HUGIN_HFOV=2
HUGIN_WIDTH=1024
HUGIN_HEIGHT=512

# options for the programs
NONA_LDR_REMAPPED_COMP=-z LZW
NONA_OPTS=
ENBLEND_OPTS= -f512x512+256+0
ENBLEND_LDR_COMP=--compression=LZW
ENBLEND_EXPOSURE_COMP=--compression=LZW
ENBLEND_HDR_COMP=
HDRMERGE_OPTS=
ENFUSE_OPTS=
EXIFTOOL_COPY_ARGS=-ImageDescription -Make -Model -Artist -WhitePoint -Copyright -GPS:all -DateTimeOriginal -CreateDate -UserComment -ColorSpace -OwnerName -SerialNumber
EXIFTOOL_INFO_ARGS='-Software=Hugin 2013.0.0.4692917e7a55' '-UserComment<$${UserComment}&\#xa;Projection: Equirectangular (2)&\#xa;FOV: 2 x 1&\#xa;Ev: -1.38' -f

# the output panorama
LDR_REMAPPED_PREFIX=auto
LDR_REMAPPED_PREFIX_SHELL=auto
HDR_STACK_REMAPPED_PREFIX=auto_hdr_
HDR_STACK_REMAPPED_PREFIX_SHELL=auto_hdr_
LDR_EXPOSURE_REMAPPED_PREFIX=auto_exposure_layers_
LDR_EXPOSURE_REMAPPED_PREFIX_SHELL=auto_exposure_layers_
PROJECT_FILE=/home/ckirst/Science/Simulation/Matlab/StemCell3D/auto.pto
PROJECT_FILE_SHELL=/home/ckirst/Science/Simulation/Matlab/StemCell3D/auto.pto
LDR_BLENDED=auto.tif
LDR_BLENDED_SHELL=auto.tif
LDR_STACKED_BLENDED=auto_fused.tif
LDR_STACKED_BLENDED_SHELL=auto_fused.tif
LDR_EXPOSURE_LAYERS_FUSED=auto_blended_fused.tif
LDR_EXPOSURE_LAYERS_FUSED_SHELL=auto_blended_fused.tif
HDR_BLENDED=auto_hdr.exr
HDR_BLENDED_SHELL=auto_hdr.exr

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
LDR_LAYERS=auto0000.tif\
auto0001.tif\
auto0002.tif\
auto0003.tif
LDR_LAYERS_SHELL=auto0000.tif\
auto0001.tif\
auto0002.tif\
auto0003.tif

# remapped images (hdr)
HDR_LAYERS=auto_hdr_0000.exr\
auto_hdr_0001.exr\
auto_hdr_0002.exr\
auto_hdr_0003.exr
HDR_LAYERS_SHELL=auto_hdr_0000.exr\
auto_hdr_0001.exr\
auto_hdr_0002.exr\
auto_hdr_0003.exr

# remapped maxval images
HDR_LAYERS_WEIGHTS=auto_hdr_0000_gray.pgm\
auto_hdr_0001_gray.pgm\
auto_hdr_0002_gray.pgm\
auto_hdr_0003_gray.pgm
HDR_LAYERS_WEIGHTS_SHELL=auto_hdr_0000_gray.pgm\
auto_hdr_0001_gray.pgm\
auto_hdr_0002_gray.pgm\
auto_hdr_0003_gray.pgm

# stacked hdr images
HDR_STACK_0=auto_stack_hdr_0000.exr
HDR_STACK_0_SHELL=auto_stack_hdr_0000.exr
HDR_STACK_0_INPUT=auto_hdr_0000.exr\
auto_hdr_0001.exr\
auto_hdr_0002.exr\
auto_hdr_0003.exr
HDR_STACK_0_INPUT_SHELL=auto_hdr_0000.exr\
auto_hdr_0001.exr\
auto_hdr_0002.exr\
auto_hdr_0003.exr
HDR_STACKS_NUMBERS=0 
HDR_STACKS=$(HDR_STACK_0) 
HDR_STACKS_SHELL=$(HDR_STACK_0_SHELL) 

# number of image sets with similar exposure
LDR_EXPOSURE_LAYER_0=auto_exposure_0000.tif
LDR_EXPOSURE_LAYER_0_SHELL=auto_exposure_0000.tif
LDR_EXPOSURE_LAYER_0_INPUT=auto_exposure_layers_0000.tif
LDR_EXPOSURE_LAYER_0_INPUT_SHELL=auto_exposure_layers_0000.tif
LDR_EXPOSURE_LAYER_0_INPUT_PTMENDER=auto0000.tif
LDR_EXPOSURE_LAYER_0_INPUT_PTMENDER_SHELL=auto0000.tif
LDR_EXPOSURE_LAYER_0_EXPOSURE=0
LDR_EXPOSURE_LAYER_1=auto_exposure_0001.tif
LDR_EXPOSURE_LAYER_1_SHELL=auto_exposure_0001.tif
LDR_EXPOSURE_LAYER_1_INPUT=auto_exposure_layers_0001.tif
LDR_EXPOSURE_LAYER_1_INPUT_SHELL=auto_exposure_layers_0001.tif
LDR_EXPOSURE_LAYER_1_INPUT_PTMENDER=auto0001.tif
LDR_EXPOSURE_LAYER_1_INPUT_PTMENDER_SHELL=auto0001.tif
LDR_EXPOSURE_LAYER_1_EXPOSURE=-3.67126
LDR_EXPOSURE_LAYER_2=auto_exposure_0002.tif
LDR_EXPOSURE_LAYER_2_SHELL=auto_exposure_0002.tif
LDR_EXPOSURE_LAYER_2_INPUT=auto_exposure_layers_0002.tif
LDR_EXPOSURE_LAYER_2_INPUT_SHELL=auto_exposure_layers_0002.tif
LDR_EXPOSURE_LAYER_2_INPUT_PTMENDER=auto0002.tif
LDR_EXPOSURE_LAYER_2_INPUT_PTMENDER_SHELL=auto0002.tif
LDR_EXPOSURE_LAYER_2_EXPOSURE=1.08515
LDR_EXPOSURE_LAYER_3=auto_exposure_0003.tif
LDR_EXPOSURE_LAYER_3_SHELL=auto_exposure_0003.tif
LDR_EXPOSURE_LAYER_3_INPUT=auto_exposure_layers_0003.tif
LDR_EXPOSURE_LAYER_3_INPUT_SHELL=auto_exposure_layers_0003.tif
LDR_EXPOSURE_LAYER_3_INPUT_PTMENDER=auto0003.tif
LDR_EXPOSURE_LAYER_3_INPUT_PTMENDER_SHELL=auto0003.tif
LDR_EXPOSURE_LAYER_3_EXPOSURE=-2.91615
LDR_EXPOSURE_LAYERS_NUMBERS=0 1 2 3 
LDR_EXPOSURE_LAYERS=$(LDR_EXPOSURE_LAYER_0) $(LDR_EXPOSURE_LAYER_1) $(LDR_EXPOSURE_LAYER_2) $(LDR_EXPOSURE_LAYER_3) 
LDR_EXPOSURE_LAYERS_SHELL=$(LDR_EXPOSURE_LAYER_0_SHELL) $(LDR_EXPOSURE_LAYER_1_SHELL) $(LDR_EXPOSURE_LAYER_2_SHELL) $(LDR_EXPOSURE_LAYER_3_SHELL) 
LDR_EXPOSURE_LAYERS_REMAPPED=auto_exposure_layers_0000.tif\
auto_exposure_layers_0001.tif\
auto_exposure_layers_0002.tif\
auto_exposure_layers_0003.tif
LDR_EXPOSURE_LAYERS_REMAPPED_SHELL=auto_exposure_layers_0000.tif\
auto_exposure_layers_0001.tif\
auto_exposure_layers_0002.tif\
auto_exposure_layers_0003.tif

# stacked ldr images
LDR_STACK_0=auto_stack_ldr_0000.tif
LDR_STACK_0_SHELL=auto_stack_ldr_0000.tif
LDR_STACK_0_INPUT=auto_exposure_layers_0000.tif\
auto_exposure_layers_0001.tif\
auto_exposure_layers_0002.tif\
auto_exposure_layers_0003.tif
LDR_STACK_0_INPUT_SHELL=auto_exposure_layers_0000.tif\
auto_exposure_layers_0001.tif\
auto_exposure_layers_0002.tif\
auto_exposure_layers_0003.tif
LDR_STACKS_NUMBERS=0 
LDR_STACKS=$(LDR_STACK_0) 
LDR_STACKS_SHELL=$(LDR_STACK_0_SHELL) 

all : startStitching $(LDR_LAYERS) 

startStitching : 
	@echo '==========================================================================='
	@echo 'Stitching panorama'
	@echo '==========================================================================='

clean : 
	@echo '==========================================================================='
	@echo 'Remove temporary files'
	@echo '==========================================================================='

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
	@echo 'Project file: /home/ckirst/Science/Simulation/Matlab/StemCell3D/auto.pto'
	@echo 'Output prefix: auto'
	@echo 'Projection: Equirectangular (2)'
	@echo 'Field of view: 2 x 1'
	@echo 'Canvas dimensions: 1024 x 512'
	@echo 'Crop area: (256,0) - (768,512)'
	@echo 'Output exposure value: -1.38'
	@echo 'Output stacks minimum overlap: 0.700'
	@echo 'Output layers maximum Ev difference: 0.50'
	@echo 'Selected outputs'
	@echo 'Normal panorama'
	@echo '* Remapped images'
	@echo '==========================================================================='
	@echo 'Input images'
	@echo '==========================================================================='
	@echo 'Number of images in project file: 4'
	@echo 'Number of active images: 4'
	@echo 'Image 0: /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0001.tif'
	@echo 'Image 0: Size 512x512, Exposure: 0.00'
	@echo 'Image 1: /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0002.tif'
	@echo 'Image 1: Size 512x512, Exposure: -3.67'
	@echo 'Image 2: /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0003.tif'
	@echo 'Image 2: Size 512x512, Exposure: 1.09'
	@echo 'Image 3: /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0004.tif'
	@echo 'Image 3: Size 512x512, Exposure: -2.92'

# Rules for ordinary TIFF_m and hdr output

auto0000.tif : /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0001.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) $(NONA_LDR_REMAPPED_COMP) -r ldr -m TIFF_m -o $(LDR_REMAPPED_PREFIX_SHELL) -i 0 $(PROJECT_FILE_SHELL)

auto_hdr_0000.exr : /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0001.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) -r hdr -m EXR_m -o $(HDR_STACK_REMAPPED_PREFIX_SHELL) -i 0 $(PROJECT_FILE_SHELL)

auto0001.tif : /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0002.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) $(NONA_LDR_REMAPPED_COMP) -r ldr -m TIFF_m -o $(LDR_REMAPPED_PREFIX_SHELL) -i 1 $(PROJECT_FILE_SHELL)

auto_hdr_0001.exr : /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0002.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) -r hdr -m EXR_m -o $(HDR_STACK_REMAPPED_PREFIX_SHELL) -i 1 $(PROJECT_FILE_SHELL)

auto0002.tif : /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0003.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) $(NONA_LDR_REMAPPED_COMP) -r ldr -m TIFF_m -o $(LDR_REMAPPED_PREFIX_SHELL) -i 2 $(PROJECT_FILE_SHELL)

auto_hdr_0002.exr : /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0003.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) -r hdr -m EXR_m -o $(HDR_STACK_REMAPPED_PREFIX_SHELL) -i 2 $(PROJECT_FILE_SHELL)

auto0003.tif : /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0004.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) $(NONA_LDR_REMAPPED_COMP) -r ldr -m TIFF_m -o $(LDR_REMAPPED_PREFIX_SHELL) -i 3 $(PROJECT_FILE_SHELL)

auto_hdr_0003.exr : /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0004.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) -r hdr -m EXR_m -o $(HDR_STACK_REMAPPED_PREFIX_SHELL) -i 3 $(PROJECT_FILE_SHELL)

# Rules for exposure layer output

auto_exposure_layers_0000.tif : /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0001.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) $(NONA_LDR_REMAPPED_COMP) -r ldr -e 0 -m TIFF_m -o $(LDR_EXPOSURE_REMAPPED_PREFIX_SHELL) -i 0 $(PROJECT_FILE_SHELL)

auto_exposure_layers_0001.tif : /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0002.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) $(NONA_LDR_REMAPPED_COMP) -r ldr -e -3.67126 -m TIFF_m -o $(LDR_EXPOSURE_REMAPPED_PREFIX_SHELL) -i 1 $(PROJECT_FILE_SHELL)

auto_exposure_layers_0002.tif : /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0003.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) $(NONA_LDR_REMAPPED_COMP) -r ldr -e 1.08515 -m TIFF_m -o $(LDR_EXPOSURE_REMAPPED_PREFIX_SHELL) -i 2 $(PROJECT_FILE_SHELL)

auto_exposure_layers_0003.tif : /home/ckirst/Science/Simulation/Matlab/StemCell3D/test0004.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) $(NONA_LDR_REMAPPED_COMP) -r ldr -e -2.91615 -m TIFF_m -o $(LDR_EXPOSURE_REMAPPED_PREFIX_SHELL) -i 3 $(PROJECT_FILE_SHELL)

# Rules for LDR and HDR stack merging, a rule for each stack

$(LDR_STACK_0) : $(LDR_STACK_0_INPUT) 
	$(ENFUSE) $(ENFUSE_OPTS) -o $(LDR_STACK_0_SHELL) -- $(LDR_STACK_0_INPUT_SHELL)
	-$(EXIFTOOL) -overwrite_original_in_place -TagsFromFile $(INPUT_IMAGE_1_SHELL) $(EXIFTOOL_COPY_ARGS) $(LDR_STACK_0_SHELL)

$(HDR_STACK_0) : $(HDR_STACK_0_INPUT) 
	$(HDRMERGE) $(HDRMERGE_OPTS) -o $(HDR_STACK_0_SHELL) -- $(HDR_STACK_0_INPUT_SHELL)

$(LDR_BLENDED) : $(LDR_LAYERS) 
	$(ENBLEND) $(ENBLEND_LDR_COMP) $(ENBLEND_OPTS) -o $(LDR_BLENDED_SHELL) -- $(LDR_LAYERS_SHELL)
	-$(EXIFTOOL) -E -overwrite_original_in_place -TagsFromFile $(INPUT_IMAGE_1_SHELL) $(EXIFTOOL_COPY_ARGS) $(EXIFTOOL_INFO_ARGS) $(LDR_BLENDED_SHELL)

$(LDR_EXPOSURE_LAYER_0) : $(LDR_EXPOSURE_LAYER_0_INPUT) 
	$(ENBLEND) $(ENBLEND_EXPOSURE_COMP) $(ENBLEND_OPTS) -o $(LDR_EXPOSURE_LAYER_0_SHELL) -- $(LDR_EXPOSURE_LAYER_0_INPUT_SHELL)
	-$(EXIFTOOL) -overwrite_original_in_place -TagsFromFile $(INPUT_IMAGE_1_SHELL) $(EXIFTOOL_COPY_ARGS) $(LDR_EXPOSURE_LAYER_0_SHELL)

$(LDR_EXPOSURE_LAYER_1) : $(LDR_EXPOSURE_LAYER_1_INPUT) 
	$(ENBLEND) $(ENBLEND_EXPOSURE_COMP) $(ENBLEND_OPTS) -o $(LDR_EXPOSURE_LAYER_1_SHELL) -- $(LDR_EXPOSURE_LAYER_1_INPUT_SHELL)
	-$(EXIFTOOL) -overwrite_original_in_place -TagsFromFile $(INPUT_IMAGE_1_SHELL) $(EXIFTOOL_COPY_ARGS) $(LDR_EXPOSURE_LAYER_1_SHELL)

$(LDR_EXPOSURE_LAYER_2) : $(LDR_EXPOSURE_LAYER_2_INPUT) 
	$(ENBLEND) $(ENBLEND_EXPOSURE_COMP) $(ENBLEND_OPTS) -o $(LDR_EXPOSURE_LAYER_2_SHELL) -- $(LDR_EXPOSURE_LAYER_2_INPUT_SHELL)
	-$(EXIFTOOL) -overwrite_original_in_place -TagsFromFile $(INPUT_IMAGE_1_SHELL) $(EXIFTOOL_COPY_ARGS) $(LDR_EXPOSURE_LAYER_2_SHELL)

$(LDR_EXPOSURE_LAYER_3) : $(LDR_EXPOSURE_LAYER_3_INPUT) 
	$(ENBLEND) $(ENBLEND_EXPOSURE_COMP) $(ENBLEND_OPTS) -o $(LDR_EXPOSURE_LAYER_3_SHELL) -- $(LDR_EXPOSURE_LAYER_3_INPUT_SHELL)
	-$(EXIFTOOL) -overwrite_original_in_place -TagsFromFile $(INPUT_IMAGE_1_SHELL) $(EXIFTOOL_COPY_ARGS) $(LDR_EXPOSURE_LAYER_3_SHELL)

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
