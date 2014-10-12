%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test ImageSource Class %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
close all
clc

initialize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ImageInfo

clc
i = ImageInfo;
i.fromData(rand(10,20))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ImageSource - base class

%% create base class with some data in memory

clear all
clear classes
close all
clc


img = ImageSource(rand(10,20))
img.setName('test');

img.printInfo

%%
clc
img.dataSize
img.dataFormat
img.dataClass
img.color
%%
clc
dat = img.data;
size(dat)


figure(1); clf; 
implot(dat)


%%
img.dataFormatPosition('Y')

%%

clc
dat = img.dataSubset('C', 1);
size(dat)

dat2 = img.dataSubset('Y', 1);
size(dat2)

figure(1); clf; 
implottiling({dat; dat2}, 'link', false)


%%
clc
figure(1); clf
imsubplot(3,1,1);
img.setName('test image');
img.setColor({'gray'});
img.plot

imsubplot(3,1,2);
img.setColor('b');
img.plot

imsubplot(3,1,3);
img.setColor('r');
img.plot


%% cells

clc
img.cellIndex(1)

%%

clc
img.cell(1)

%%
clc
img.data


%% roi

roi = ROIRectangle([5,5], [10,8]);

dd = img.dataExtract(roi)

size(dd)


%%%%%%%%%%%%%%%%%%%%%%%%%
%% ImageSourceFile


clear all
clear classes
close all
clc

initialize
bfinitialize


%% accessing a single file
clc
is = ImageSourceFile('./Test/Images/hESCells_DAPI.tif');
is.printInfo

%%

is.setCache(false);

size(is.idata)

img = is.data;

figure(1)
is.setColor('r');
is.plot

% no caching 
size(is.idata)


%%

is.setCache(true);

size(is.idata)

img = is.data;

figure(1)
is.setColor('r');
is.plot

% now dat is cached
size(is.idata)


%%
is.clearCache


%% raw vs data

clear all
clear classes
close all
clc

initialize

%%
clc
is = ImageSourceFile('./Test/Images/hESCells_DAPI.tif');
is.setRawDataFormat('Xy');
is.setCache(false);

is

%%

figure(1)
imsubplot(2,1,1)
is.setColor('r');
is.setRawDataFormat('XY');
is.plot

imsubplot(2,1,2)
is.setColor('b');
is.setRawDataFormat('Xy');
is.plot

