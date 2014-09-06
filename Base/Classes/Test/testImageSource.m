%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test ImageSource Class %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
close all
clc

%% Basic Functionality

%% create base class with some data in memory

img = ImageSource(rand(10,20))
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
img.setColor('red')
img.plot


figure(3)
img.setColor('blue')
img.plot


%%%%%%%%%%%%%%%%%%%%%%%%%
%% Labeled data access


clear all
clear classes
close all
clc


%%

ImageFormatIndex('c', 1)

%%
d = ImageFormatIndex('c', 1, 't', 7)

{d.icoordinate}
{d.iindex}

%%
clc

ifl = ImageFormatLabel('dapi', {'c', 2, 't', 2:4})

ifi = ifl.ilabel('dapi')
ifi(2)


ifi.indices('pqct')


%%

img = ImageSource(rand(10,20))

ifl = ImageFormatLabel('dapi', {'p', 5})
img.setLabel(ifl)

img.label

%%
img.extract('dapi') - img.idata(5,:)


%% 
