%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test ImageSourceTiled Class %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
close all
clc

initialize
bfinitialize


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ImageSourceTiled 

%% prepare source

clc
texpr = tagexpr(dirr('./Test/Images/hESCells_Tiling/*.tif'), 'tagnames', {'tile'})
is = ImageSourceTagged(texpr);
is.print

is.itagranges

%% reudce to first 4 images

is.setTagRange('tile',  {37, 38, 33, 34})
is.print
is.cellsize

%% test id data access
clc
d = is.data('tile', 38);

d2 = is.data(2);

figure(1);
implottiling({d; d2});


%%

clc
c = is.celldata;
size(c)

figure(1);
implottiling(c)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ImageSourceTiled - more testing

clear all
clear classes
close all
clc

initialize
bfinitialize


%%
clc

texpr = tagexpr(dirr('./Test/Images/hESCells_Tiling/*.tif'), 'tagnames', {'tile'})
is = ImageSourceTagged(texpr);
is.tagrangesize

%%

is.setTagRange('tile',  {45, 46, 47, 48,  41, 42, 43, 44,  37, 38, 39, 40,  33, 34, 35, 36})
ist = ImageSourceTiled(is, 'tileshape', [4,4])


%% get Tiles

tl = ist.tiles

figure(1);
implottiling(tl, 'tiling', [4,4])


%%
ts = ist.getTileSizes
ts{1}

%%

ist.ntiles

