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
%% ImageSourceFiles - init

clc
is = ImageSourceFiles
is.printInfo


%%

is.initializeRawCellFromFileExpression('./Test/Images/hESCells_Stack/W<W,1>F<F,3>T0001Z<Z,2>C1.tif')
is.printInfo


%%

is.initializeRawDataFromFile(is.fileNameFromRawRange(1))
is.printInfo


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ImageSourceFiles

clear all
clear classes
close all
clc

initialize
bfinitialize

obj = ImageSourceFiles;


%%

clc
is = ImageSourceFiles('./Test/Images/hESCells_Stack/W<W,1>F<F,3>T0001Z<Z,2>C<C,1>.tif');
is.printInfo

%%

is.fileTagRange

%%

is.rawCellDataFormat
is.rawCellDataSize

%%
clc
is.fileName('W', 1)

%%
clc
is.fileName('Z', 10)

%%
clc
is.fileName('z', 10)

%%

is.fileName('z', 1:5, 'W', 1)

%%

is.fileName

%%

is.fileName(1,1,6,1)


%% 

is.setKey('Z', num2cell('abcdefghiklmnop'))
is.key
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
is.fileTagRange

is.rangeFromFileTagRange

%%

is.setRange('F', 2:5)

is.fileTagRange


%%
clc
is.rawRangeFromFileTagRange


%%
clc
is.rangeFromFileTagRange('F', 34)

%%

clc
is.fileName(1:3)

%%

is.range

is.cellDataFormat
is.cellDataSize

%%

is.fileName('F', 1:5, 'W', 1)
is.fileNameFromRawRange

%%

d = is.data(1);
size(d)
is.dataSize

%%
figure(1); clf;
implot(d)

%%
clc
cdat = is.cell;

size(cdat)
is.cellSize


%%

is.setRawCellDataCaching(true);

cdat = is.cell

figure(1); clf
implottiling(cdat')

%%

is.setRawCellDataCaching(false);

%%

is.initializeDataIntensity;
is.maxDataIntensity
is.minDataIntensity


%%
clc
is.resetRange
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



