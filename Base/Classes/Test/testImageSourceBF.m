%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test ImageSourceBF Class %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
close all
clc

initialize
bfinitialize

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
size(is.data)

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


%% multi channel data


clear all
clear classes
close all
clc

initialize


is = ImageSourceBF('./Test/Images/hESCells_Colony.zvi')


%%
clc
[ts,tf]= is.tileSizeAndFormat


%% reshape
is.setReshape('S', 'UV', ts);
is.setCellFormat(tf);
is.setRange('C', 1)

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



