%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test ImageSourceTagged Class %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
close all
clc

initialize
bfinitialize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ImageSourceTagged - base class

clc
is = ImageSourceTagged('./Test/Images/hESCells_Stack/W<well,1>F<field,3>T0001Z<z,2>C1.tif');
is.print


%%

is.filename('z', 50)
is.command('z', 50)
is.infocommand('z', 1, 'field', 127, 'well', 1)

%%

is.tag(6)
is.tagrangesize

%%

clc
is.datasize
is.cellsize

is.dataformat

is.color

%% getInfo

clc
is.getInfo


%%

is.print

%is.initialze();


%% 

clc
is.dataformat

[i, p] = is.rawtagids

[i,p] = is.rawtagidsinternal
[i,p] = is.rawtagidsexternal


%%
is.celltagranges()


%%
clc
is.rawtagsranges
is.rawtagrangesinternal
is.rawtagrangesexternal
is.celltagranges

%%
clc
tgrs.z = num2cell([1,2,3,5,6]);

is.rawtagsranges(tgrs)
is.parseTagRanges(tgrs)

%% data
clc
d = is.data('z', {6,7,8,9});

size(d)
imformat(d)

figure(1); clf
implottiling(mat2gray(d))



%% celldata

cd = is.celldata('z', {1,2,3})

size(cd)


%%
is.celltagnames
is.datatagnames



%% Test some Tagrange manipulations

ir = is.itagranges

ir.tag1 = {1,4};
ir.p = {1,2,3,4};
ir.tag2 = {'a', 'b', 'd'};

is.itagranges   = ir;
is.itagformat   = is.tagformatFromTagRanges
is.itaginternal = [0, 0, 1,0]


%%
is.print


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% different read commands

clear all
clear classes
close all
clc

initialize
bfinitialize

%%

clc
is = ImageSourceTagged('filename', './Test/Images/hESCells_Stack/W1F<field,3>T0001Z<z,2>C1.tif', 'readcommand', 'imreadBF(''<file>'', ''p'', <p>)', ...
   'tagranges', struct('p', {num2cell(1:50)}));
is.print

%%
is.tagrangesize


%%

d =is.data('z', {1,5,10,15});

figure(1);
implottiling(d)

%%
is.datasize('z', {1,5,10,15})
is.dataformat
size(d)

%%

is.dataformat('z', {3,4}, 'p', 1)

%%

dat = is.data('z', 1);
img = imreadBF('./Test/Images/hESCells_Stack/W1F127T0001Z01C1.tif', 'p', 1:50);

figure(1);
implottiling({img; dat} )


