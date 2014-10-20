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

is.metaDataNames

is.metaData('XResolution')


%%
clc
is.rangeNames
is.nCells


%%
clc
is.rawRange

%%

is.setRange('S', 1);
is.printInfo


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
is.rawCellDataCache


%%
clc
is.setRawCaching(true);
is.clearCache

var2char({size(is.cellDataCache), size(is.rawCellDataCache)})

figure(1)
is.plot

% should notgive asecond read
is.data(1);

% now dat is cached on raw files
var2char({size(is.cellDataCache), size(is.rawCellDataCache)})


%%
clc
is.setRawCaching(false);
is.clearCache

var2char({size(is.cellDataCache), size(is.rawCellDataCache)})

figure(1)
is.plot

%should give a second read
is.data(1);

% dat is not cached on raw files
var2char({size(is.cellDataCache), size(is.rawCellDataCache)})


%%
is.clearCache;

%% raw vs data

clear all
clear classes
close all
clc

initialize

%%
clc
is = ImageSourceBF('./Test/Images/hESCells_DAPI.tif');

%%
clc


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

is.key
is.range('C', 'GFP')


%%

clc
d =is.data('C', 'DAPI', 'S', 1);

size(d)


%%
clc
[ts,tf]= is.tileSizeAndFormat


%% reshape
is.setReshape('S', 'UV', ts);
is.setCellFormat(tf);
is.setRange('C', 1);

is.printInfo

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
is.setCaching(false);

is.initializeFromTiling;
is.setRange('C', 1);
is.printInfo


figure(1); clf
cd = is.cell;
implottiling(cd)

%%
clc
clf
is.cellPreview('overlap', 140)


%%

is.cellFormat


%%
is.rangeFromCellIndex(1,2)

%%
clc
is.dataSize('C', 1)
is.dataSize
is.dataSizeFromRawSize
is.fullDataSize

%%
clc
is.cellSize('U', 1)
is.cellSize
is.cellSizeFromRawSize
is.fullCellSize

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

bkg = zeros(is.dataSize);
size(bkg)
flt = 0.5 *  ones(is.dataSize);

is.setBackgroundAndFlatFieldCorrection(bkg, flt);

is.dataCorrect


%%

is.setDataCorrect(true);
ds = is.correctData(ones(is.dataSize));
size(ds)
max(ds(:))

is.setDataCorrect(false);
ds = is.correctData(ones(is.dataSize));
size(ds)
max(ds(:))


%%
is.clearCache
is.setCaching(false)

is.setDataCorrect(true);
d = is.data(1);

is.setDataCorrect(false);
d2 = is.data(1);

figure(2); clf
implottiling({d; d2})

d(1) - d2(1)


%% arbitrary backgroudn correction


is.setDataCorrectFunction(@(x) x*2);

is.setDataCorrect(true);
d = is.data(1);

is.setDataCorrect(false);
d2 = is.data(1);

figure(2); clf
implottiling({d; d2})



%% keys

clear all
clear classes
close all
clc

initialize

obj = ImageSourceBF;

%%
is = ImageSourceBF('./Test/Images/hESCells_Colony.zvi');
is.setCaching(true);

is.initializeFromTiling;
is.setRange('C', 1);
is.printInfo

is.key

%%
is.dataFormat

is.addKey('X', cellfunc(@num2str, num2cell(0:1377)));
is.key

%%
is.range('C', 'GFP', 'X', '0')

%%


is.rawRange('C', 'GFP', 'X', '0')


%%
clc
is.rawRangeFromRawRange('c', 1)




