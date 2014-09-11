%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test ImageSource Class %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
close all
clc

initialize

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
size(img)
format(img)
class(img)
color(img)

img.sdims

dat = data(img);
size(dat)

%%
clc
figure(1);
img.setName('test image');
img.plot

figure(2)
img.setColor('gray')
img.plot


figure(3)
img.setColor('blue')
img.plot




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
is.initialize;

is

%%

is.print


%%

is.icache = false;

size(is.idata)

img = is.data;

figure(1)
is.setColor('r');
is.plot

% no caching 
size(is.idata)


%%

is.icache = true;

size(is.idata)

img = is.data;

figure(1)
is.setColor('r');
is.plot

% now dat is cached
size(is.idata)


















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


