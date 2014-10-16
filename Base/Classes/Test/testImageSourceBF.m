%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test ImageSourceBF Class %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
close all
clc

initialize
bfinitialize

obj = ImageSource

%%

ii = ImageInfo('./Test/Images/hESCells_DAPI.tif')

%% accessing a single file
clc
is = ImageSourceBF('./Test/Images/hESCells_DAPI.tif');
is.printInfo

%%
md = is.metaData

is.metaData('XResolution')

%%
clc
is.rangeNames
is.nCells


%%
clc
is.rawRange

%%

is.setRange('S', 1)


%%

is.setCaching(false);

%% data access
clear all
clear classes
close all
clc

initialize
bfinitialize

%%
clc
is = ImageSourceBF('./Test/Images/hESCells_DAPI.tif');
is.printInfo

%%

d = is.data;
size(d)

%%
clc
figure(1)
is.plot

%% no caching 

is.cellDataCache


%%
clc
is.setCaching(true);

size(is.cellDataCache)

img = is.data;

figure(1)
is.plot

% now dat is cached
size(is.cellDataCache)


%%
is.clearCellDataCache;

%% raw vs data

clear all
clear classes
close all
clc

initialize

%%
clc
is = ImageSourceBF('./Test/Images/hESCells_DAPI.tif');
is.setCaching(false);

is

%%
clc
is.setCaching(false);

figure(1); clf
imsubplot(2,1,1)
is.setDataFormat('XY');
d = is.data;
implot(d)

imsubplot(2,1,2)
is.setDataFormat('yX')
d=is.data;
implot(d)



is.rawDataFormat
is.dataFormat


%%

figure(1); clf
imsubplot(2,1,1)
is.setDataFormat('XY');
is.plot

imsubplot(2,1,2)
is.setDataFormat('Xy');
is.plot



%%
clc

is.printInfo


%% multi channel / tiled data


clear all
clear classes
close all
clc

initialize
bfinitialize

%%
is = ImageSourceBF('./Test/Images/hESCells_Colony.zvi')

is.printInfo

%% channel keys

is.rangeKey
is.rangeToIndexRange('C', 'GFP')


%%

clc
d =is.data('C', 'DAPI', 'S', 1);

size(d)

%%
is.rangeKey


%%
clc
[ts,tf]= is.tileSizeAndFormat


%% reshape
is.setReshape('S', 'UV', ts);
is.setCellFormat(tf);
is.setRange('C', 1)

%%

is.dataSize
is.cellSize

%%
clc
figure(1); clf
cd = is.cell;
implottiling(cd)


%% 

clear all
clear classes
close all
clc

initialize

%%
is = ImageSourceBF('./Test/Images/hESCells_Colony.zvi');
is.initializeFromTiling;
is.setRange('C', 1);
is.printInfo


figure(1); clf
cd = is.cell;
implottiling(cd)

%%
clc
is.previewCellAlignment('overlap', 140)



%%

is.cellIndexToCoordinate(2,2)

%%
clc
is.dataSize('C', 1)

%%
is.dataSizeFromRaw
is.cellSizeFromRaw

%%
clc
is.cellSize
is.cellFormat
is.cellSize('V', 1)

%%

d = is.data(3);
size(d)

figure(2); clf
implot(d)


%% intensities

clc
is.initializeDataIntensity('C', 1, 'U', 1);
is.minDataIntensity
is.maxDataIntensity


%% background correction

clc
is.dataSize

bkg = zeros(1378, 1024);
size(bkg)
flt = 0.5 * ones(1378, 1024);

is.setBackgroundAndFlatFieldCorrection(bkg, flt);

d = is.data(1);

is.setDataCorrect(false);

d2 = is.data(1);

figure(2); clf
implottiling({d; d2})


%% arbitrary backgroudn correction


is.setDataCorrectFunction(@(x) x*2);

is.setDataCorrect(true);
d = is.data(1);

is.setDataCorrect(false);
d2 = is.data(1);

figure(2); clf
implottiling({d; d2})












