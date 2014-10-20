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

obj = ImageSource


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






%% reformatting and ranges


clear all
clear classes
close all
clc

img = ImageSource(rand(10,20));
img.printInfo

img.setName('test');
img.setRange('Y', [1,2,5,6]);

img.printInfo

%%
clc
img.setReshape('Y', 'UV', [4,5]);
img.printInfo

%%

img.dataSize
img.cellSize

img.cellDataSize
img.cellDataFormat

%%

img.rangeFromVarargin(1)

img.rawRangeFromRawVarargin(1)

%%

img.rawCellDataSize




%%
img.range
img.rawRange



%%

d = img.data(1);
size(d)


%% 
img.addRange('V', 2);

d = img.data(1);
size(d)


img.rawRange


