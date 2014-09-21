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

clc

img = ImageSource(rand(10,20))
img.setName('test');
img.info

img.print

img.info

%%
clc
img.datasize
img.dataformat
img.dataclass
img.color

img.nsdatadims

dat = data(img);
size(dat)


img.dataformatpos('p')

%%

clc
dat = img.subdata('c', 1);
size(dat)

dat = img.subdata('p', 1);
size(dat)


%%
clc
figure(1); clf
imsubplot(3,1,1);
img.setName('test image');
img.setColor({'gray'})
img.plot

imsubplot(3,1,2);
img.setColor('b');
img.plot

imsubplot(3,1,3);
img.setColor('r');
img.plot


%%

img.data


%% roi

roi = ROIRectangle([5,5], [10,8]);

dd = img.extractdata(roi)

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

is.info

%%

is.print


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













%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data acess via Keywords: todo


clear all
clear classes
close all
clc


%%

ImageDataIndex('c', 1)

%%
d = ImageDataIndex('c', 1, 't', 7)

{d.icoordinate}
{d.iindex}

%%
clc

idl = ImageDataLabel('dapi', {'c', 2, 't', 2:4})

idi = idl.ilabel('dapi')
idi(2)

idi.indices('pqct')

idi.infoString


%%

img = ImageSource(rand(10,20))

idl = ImageDataLabel('dapi', {'p', 5}, 'hello', {'p', 4:6, 'q', 7:9})

img.setLabel(idl)

img.label.infoString

img.print

%%
img.extract('dapi') - img.idata(5,:)

img.extract('hello') - img.idata(4:6, 7:9)


