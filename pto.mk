
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
HUGIN_HFOV=360
HUGIN_WIDTH=3000
HUGIN_HEIGHT=1500

# options for the programs
NONA_LDR_REMAPPED_COMP=-z LZW
NONA_OPTS=
ENBLEND_OPTS= -w -f3000x1500
ENBLEND_LDR_COMP=--compression=LZW
ENBLEND_EXPOSURE_COMP=--compression=LZW
ENBLEND_HDR_COMP=
HDRMERGE_OPTS=
ENFUSE_OPTS= -w
EXIFTOOL_COPY_ARGS=-ImageDescription -Make -Model -Artist -WhitePoint -Copyright -GPS:all -DateTimeOriginal -CreateDate -UserComment -ColorSpace -OwnerName -SerialNumber
EXIFTOOL_INFO_ARGS='-Software=Hugin 2013.0.0.4692917e7a55' '-UserComment<$${UserComment}&\#xa;Projection: Equirectangular (2)&\#xa;FOV: 360 x 180&\#xa;Ev: -1.55' -f

# the output panorama
LDR_REMAPPED_PREFIX=prefix
LDR_REMAPPED_PREFIX_SHELL=prefix
HDR_STACK_REMAPPED_PREFIX=prefix_hdr_
HDR_STACK_REMAPPED_PREFIX_SHELL=prefix_hdr_
LDR_EXPOSURE_REMAPPED_PREFIX=prefix_exposure_layers_
LDR_EXPOSURE_REMAPPED_PREFIX_SHELL=prefix_exposure_layers_
PROJECT_FILE=auto.pto
PROJECT_FILE_SHELL=auto.pto
LDR_BLENDED=prefix.tif
LDR_BLENDED_SHELL=prefix.tif
LDR_STACKED_BLENDED=prefix_fused.tif
LDR_STACKED_BLENDED_SHELL=prefix_fused.tif
LDR_EXPOSURE_LAYERS_FUSED=prefix_blended_fused.tif
LDR_EXPOSURE_LAYERS_FUSED_SHELL=prefix_blended_fused.tif
HDR_BLENDED=prefix_hdr.exr
HDR_BLENDED_SHELL=prefix_hdr.exr

# first input image
INPUT_IMAGE_1=test0001.tif
INPUT_IMAGE_1_SHELL=test0001.tif

# all input images
INPUT_IMAGES=test0001.tif\
test0002.tif\
test0003.tif\
test0004.tif
INPUT_IMAGES_SHELL=test0001.tif\
test0002.tif\
test0003.tif\
test0004.tif

# remapped images
LDR_LAYERS=prefix0000.tif\
prefix0001.tif\
prefix0002.tif\
prefix0003.tif
LDR_LAYERS_SHELL=prefix0000.tif\
prefix0001.tif\
prefix0002.tif\
prefix0003.tif

# remapped images (hdr)
HDR_LAYERS=prefix_hdr_0000.exr\
prefix_hdr_0001.exr\
prefix_hdr_0002.exr\
prefix_hdr_0003.exr
HDR_LAYERS_SHELL=prefix_hdr_0000.exr\
prefix_hdr_0001.exr\
prefix_hdr_0002.exr\
prefix_hdr_0003.exr

# remapped maxval images
HDR_LAYERS_WEIGHTS=prefix_hdr_0000_gray.pgm\
prefix_hdr_0001_gray.pgm\
prefix_hdr_0002_gray.pgm\
prefix_hdr_0003_gray.pgm
HDR_LAYERS_WEIGHTS_SHELL=prefix_hdr_0000_gray.pgm\
prefix_hdr_0001_gray.pgm\
prefix_hdr_0002_gray.pgm\
prefix_hdr_0003_gray.pgm

# stacked hdr images
HDR_STACK_0=prefix_stack_hdr_0000.exr
HDR_STACK_0_SHELL=prefix_stack_hdr_0000.exr
HDR_STACK_0_INPUT=prefix_hdr_0000.exr\
prefix_hdr_0001.exr\
prefix_hdr_0002.exr\
prefix_hdr_0003.exr
HDR_STACK_0_INPUT_SHELL=prefix_hdr_0000.exr\
prefix_hdr_0001.exr\
prefix_hdr_0002.exr\
prefix_hdr_0003.exr
HDR_STACKS_NUMBERS=0 
HDR_STACKS=$(HDR_STACK_0) 
HDR_STACKS_SHELL=$(HDR_STACK_0_SHELL) 

# number of image sets with similar exposure
LDR_EXPOSURE_LAYER_0=prefix_exposure_0000.tif
LDR_EXPOSURE_LAYER_0_SHELL=prefix_exposure_0000.tif
LDR_EXPOSURE_LAYER_0_INPUT=prefix_exposure_layers_0000.tif
LDR_EXPOSURE_LAYER_0_INPUT_SHELL=prefix_exposure_layers_0000.tif
LDR_EXPOSURE_LAYER_0_INPUT_PTMENDER=prefix0000.tif
LDR_EXPOSURE_LAYER_0_INPUT_PTMENDER_SHELL=prefix0000.tif
LDR_EXPOSURE_LAYER_0_EXPOSURE=0
LDR_EXPOSURE_LAYER_1=prefix_exposure_0001.tif
LDR_EXPOSURE_LAYER_1_SHELL=prefix_exposure_0001.tif
LDR_EXPOSURE_LAYER_1_INPUT=prefix_exposure_layers_0001.tif
LDR_EXPOSURE_LAYER_1_INPUT_SHELL=prefix_exposure_layers_0001.tif
LDR_EXPOSURE_LAYER_1_INPUT_PTMENDER=prefix0001.tif
LDR_EXPOSURE_LAYER_1_INPUT_PTMENDER_SHELL=prefix0001.tif
LDR_EXPOSURE_LAYER_1_EXPOSURE=-3.97454
LDR_EXPOSURE_LAYER_2=prefix_exposure_0002.tif
LDR_EXPOSURE_LAYER_2_SHELL=prefix_exposure_0002.tif
LDR_EXPOSURE_LAYER_2_INPUT=prefix_exposure_layers_0002.tif
LDR_EXPOSURE_LAYER_2_INPUT_SHELL=prefix_exposure_layers_0002.tif
LDR_EXPOSURE_LAYER_2_INPUT_PTMENDER=prefix0002.tif
LDR_EXPOSURE_LAYER_2_INPUT_PTMENDER_SHELL=prefix0002.tif
LDR_EXPOSURE_LAYER_2_EXPOSURE=1.07964
LDR_EXPOSURE_LAYER_3=prefix_exposure_0003.tif
LDR_EXPOSURE_LAYER_3_SHELL=prefix_exposure_0003.tif
LDR_EXPOSURE_LAYER_3_INPUT=prefix_exposure_layers_0003.tif
LDR_EXPOSURE_LAYER_3_INPUT_SHELL=prefix_exposure_layers_0003.tif
LDR_EXPOSURE_LAYER_3_INPUT_PTMENDER=prefix0003.tif
LDR_EXPOSURE_LAYER_3_INPUT_PTMENDER_SHELL=prefix0003.tif
LDR_EXPOSURE_LAYER_3_EXPOSURE=-3.30112
LDR_EXPOSURE_LAYERS_NUMBERS=0 1 2 3 
LDR_EXPOSURE_LAYERS=$(LDR_EXPOSURE_LAYER_0) $(LDR_EXPOSURE_LAYER_1) $(LDR_EXPOSURE_LAYER_2) $(LDR_EXPOSURE_LAYER_3) 
LDR_EXPOSURE_LAYERS_SHELL=$(LDR_EXPOSURE_LAYER_0_SHELL) $(LDR_EXPOSURE_LAYER_1_SHELL) $(LDR_EXPOSURE_LAYER_2_SHELL) $(LDR_EXPOSURE_LAYER_3_SHELL) 
LDR_EXPOSURE_LAYERS_REMAPPED=prefix_exposure_layers_0000.tif\
prefix_exposure_layers_0001.tif\
prefix_exposure_layers_0002.tif\
prefix_exposure_layers_0003.tif
LDR_EXPOSURE_LAYERS_REMAPPED_SHELL=prefix_exposure_layers_0000.tif\
prefix_exposure_layers_0001.tif\
prefix_exposure_layers_0002.tif\
prefix_exposure_layers_0003.tif

# stacked ldr images
LDR_STACK_0=prefix_stack_ldr_0000.tif
LDR_STACK_0_SHELL=prefix_stack_ldr_0000.tif
LDR_STACK_0_INPUT=prefix_exposure_layers_0000.tif\
prefix_exposure_layers_0001.tif\
prefix_exposure_layers_0002.tif\
prefix_exposure_layers_0003.tif
LDR_STACK_0_INPUT_SHELL=prefix_exposure_layers_0000.tif\
prefix_exposure_layers_0001.tif\
prefix_exposure_layers_0002.tif\
prefix_exposure_layers_0003.tif
LDR_STACKS_NUMBERS=0 
LDR_STACKS=$(LDR_STACK_0) 
LDR_STACKS_SHELL=$(LDR_STACK_0_SHELL) 
DO_LDR_BLENDED=1

all : startStitching $(LDR_BLENDED) 

startStitching : 
	@echo '==========================================================================='
	@echo 'Stitching panorama'
	@echo '==========================================================================='

clean : 
	@echo '==========================================================================='
	@echo 'Remove temporary files'
	@echo '==========================================================================='
	-$(RM) $(LDR_LAYERS_SHELL) 

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
	@echo 'Project file: auto.pto'
	@echo 'Output prefix: prefix'
	@echo 'Projection: Equirectangular (2)'
	@echo 'Field of view: 360 x 180'
	@echo 'Canvas dimensions: 3000 x 1500'
	@echo 'Crop area: (0,0) - (3000,1500)'
	@echo 'Output exposure value: -1.55'
	@echo 'Output stacks minimum overlap: 0.700'
	@echo 'Output layers maximum Ev difference: 0.50'
	@echo 'Selected outputs'
	@echo 'Normal panorama'
	@echo '* Blended panorama'
	@echo '==========================================================================='
	@echo 'Input images'
	@echo '==========================================================================='
	@echo 'Number of images in project file: 4'
	@echo 'Number of active images: 4'
	@echo 'Image 0: test0001.tif'
	@echo 'Image 0: Size 512x512, Exposure: 0.00'
	@echo 'Image 1: test0002.tif'
	@echo 'Image 1: Size 512x512, Exposure: -3.97'
	@echo 'Image 2: test0003.tif'
	@echo 'Image 2: Size 512x512, Exposure: 1.08'
	@echo 'Image 3: test0004.tif'
	@echo 'Image 3: Size 512x512, Exposure: -3.30'

# Rules for ordinary TIFF_m and hdr output

prefix0000.tif : test0001.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) $(NONA_LDR_REMAPPED_COMP) -r ldr -m TIFF_m -o $(LDR_REMAPPED_PREFIX_SHELL) -i 0 $(PROJECT_FILE_SHELL)

prefix_hdr_0000.exr : test0001.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) -r hdr -m EXR_m -o $(HDR_STACK_REMAPPED_PREFIX_SHELL) -i 0 $(PROJECT_FILE_SHELL)

prefix0001.tif : test0002.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) $(NONA_LDR_REMAPPED_COMP) -r ldr -m TIFF_m -o $(LDR_REMAPPED_PREFIX_SHELL) -i 1 $(PROJECT_FILE_SHELL)

prefix_hdr_0001.exr : test0002.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) -r hdr -m EXR_m -o $(HDR_STACK_REMAPPED_PREFIX_SHELL) -i 1 $(PROJECT_FILE_SHELL)

prefix0002.tif : test0003.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) $(NONA_LDR_REMAPPED_COMP) -r ldr -m TIFF_m -o $(LDR_REMAPPED_PREFIX_SHELL) -i 2 $(PROJECT_FILE_SHELL)

prefix_hdr_0002.exr : test0003.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) -r hdr -m EXR_m -o $(HDR_STACK_REMAPPED_PREFIX_SHELL) -i 2 $(PROJECT_FILE_SHELL)

prefix0003.tif : test0004.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) $(NONA_LDR_REMAPPED_COMP) -r ldr -m TIFF_m -o $(LDR_REMAPPED_PREFIX_SHELL) -i 3 $(PROJECT_FILE_SHELL)

prefix_hdr_0003.exr : test0004.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) -r hdr -m EXR_m -o $(HDR_STACK_REMAPPED_PREFIX_SHELL) -i 3 $(PROJECT_FILE_SHELL)

# Rules for exposure layer output

prefix_exposure_layers_0000.tif : test0001.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) $(NONA_LDR_REMAPPED_COMP) -r ldr -e 0 -m TIFF_m -o $(LDR_EXPOSURE_REMAPPED_PREFIX_SHELL) -i 0 $(PROJECT_FILE_SHELL)

prefix_exposure_layers_0001.tif : test0002.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) $(NONA_LDR_REMAPPED_COMP) -r ldr -e -3.97454 -m TIFF_m -o $(LDR_EXPOSURE_REMAPPED_PREFIX_SHELL) -i 1 $(PROJECT_FILE_SHELL)

prefix_exposure_layers_0002.tif : test0003.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) $(NONA_LDR_REMAPPED_COMP) -r ldr -e 1.07964 -m TIFF_m -o $(LDR_EXPOSURE_REMAPPED_PREFIX_SHELL) -i 2 $(PROJECT_FILE_SHELL)

prefix_exposure_layers_0003.tif : test0004.tif $(PROJECT_FILE) 
	$(NONA) $(NONA_OPTS) $(NONA_LDR_REMAPPED_COMP) -r ldr -e -3.30112 -m TIFF_m -o $(LDR_EXPOSURE_REMAPPED_PREFIX_SHELL) -i 3 $(PROJECT_FILE_SHELL)

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
