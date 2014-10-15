%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test ImageSourceTagged Class %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
close all
clc

initialize
bfinitialize

obj = ImageSource;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ImageSourceFiles - base class

clc
is = ImageSourceFiles
is.printInfo


%%

is.initializeRawCellFromFileExpression('./Test/Images/hESCells_Stack/W<W,1>F<F,3>T0001Z<Z,2>C1.tif')
is.printInfo

%%
clc
is.rawFileName('W', 1)

%%
clc
is.rawFileName('Z', 10)

%%

is.initializeRawDataFromFile(is.rawFileName(1))
is.printInfo



%% Test

clear all
clear classes
close all
clc

initialize
bfinitialize

obj = ImageSource;


%%

clc
is = ImageSourceFiles('./Test/Images/hESCells_Stack/W<W,1>F<F,3>T0001Z<Z,2>C<C,1>.tif');
is.printInfo

%%

is.fileName('z', 1:5, 'W', 1)

%%

is.fileName

%%

is.fileName(1,1,2,1)


%% 

is.setRangeKey('Z', num2cell('abcdefghiklmnop'))
is.rangeKey
is.fileName('Z', 'f', 'W', 1)


%% Tilings


clear all
clear classes
close all
clc

initialize
bfinitialize

obj = ImageSource;

is = ImageSourceFiles('./Test/Images/hESCells_Tiling/W<W,1>F<F,3>T0001Z<Z,2>C<C,1>.tif');
is.printInfo

%%

d = is.data(1);
figure(1); clf;
implot(d)

%%
clc
is.cell


%%

cdat = is.cell(1:4)

figure(1); clf
implottiling(cdat')

%%

is.initializeDataIntensity;
is.maxDataIntensity
is.minDataIntensity


%%
clc
is.setReshape('F', 'UV', [4,4]);
is.reformatCellFormat('Uv');
is.printInfo

%%
cd = is.cell;
size(cd)

figure(1)
implottiling(cd)


%%
is.fileName

%%
clc
is.dataSize('X', 10:15)

%%

is.cellSize('U', 1)


%%

is.rawCellSize

%%

is.rawDataFormat
is.rawCellFormat
is.rawCellDataFormat

%%

is.range





%%
clc
p = struct('x', 1, 'y', 3:10);
p.x = {'a', 'b', 'c'};

imfrmtRangeToIndexRange(p, p)



