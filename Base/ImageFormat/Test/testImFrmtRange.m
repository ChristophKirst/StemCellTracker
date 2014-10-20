%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Image Formats    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all
initialize


%% Ranges


imfrmtRangeFromVarargin('X', 10, 'y', 1:4, 'U', 10, 'x', 5)

imfrmtRangeFromFormatAndVarargin('XYZ', 'X', 10, 'y', 1:4, 'U', 10, 'x', 5)

imfrmtRangeFromSizeFormatAndVarargin([10,20,30], 'XYZ', 'X', 10, 'y', 1:4, 'U', 10, 'x', 5)

imfrmtRangeFromRangeAndVarargin(struct('X', 1:10, 'y', 3:6), 'X', 4, 'Y', 1:2, 'u', 1, 'x', 3)


%%

imfrmtReformatRange([10,20,30], 'XYZ', 'yZ',  struct('X', 1:10, 'Y', 3:6))


%%
refRange = struct('X', {num2cell(0:10)}, 'Y', {num2cell(3:6)}, 'C', {{'a', 'b', 'c'}})


imfrmtRangeFromIndexRange(refRange, struct('X', 1:2, 'y', 2, 'C', 2, 'U', 4))

%%
clc
imfrmtRangeToIndexRange(refRange, struct('X',{num2cell(1:2)}, 'Y', {{5}}, 'C', {{'a', 'b'}}, 'U', 3))


%%

imfrmtRangeFromIndex([10,20,30], 'XYZ', 3, 2:4, 2)

%%


imfrmtIndexRangeFromIndexRange( struct('X', 0:2, 'y', 2:6, 'C', 2, 'U', 4),  struct('X', 1:2, 'y', 2, 'V', 6))



%%
clc
rg = struct('X', [3:5], 'Y', 4:5);

idx = imfrmtIndexToFullIndex([10,20], 'XY', rg, [1:6])

%%
imfrmtFullIndexToIndex([10,20], 'XY', rg, idx)


%% reshaping ranges
clc

rawSi = [10,20];
rawRg = struct('X', 3:6, 'Y', 5:12);
rawFrmt= 'XY';

frmt = 'XYZ';

refrom = {'Y'};
reto = {'YZ'};
resi = {[4,5]};

rg = imfrmtReshapeRange(rawSi, rawFrmt, frmt, refrom, reto, resi, rawRg)
si  = imfrmtReshapeSize(rawSi, rawFrmt, frmt, refrom, reto, resi);

siRange = imfrmtRangeSize(si, frmt, rg)

%%
clc
[rgsi, rshi] = imfrmtReshapeInverseRange(si, frmt, rawFrmt, refrom, reto, resi, rg);
rgsi
var2char(rshi)



%%

% test on data

rawDat = rand(10,20);

dat = imfrmtReshape(rawDat, rawFrmt, frmt, refrom, reto, resi);
size(dat)
si



%%

datRange =imfrmtDataSubset(dat, 'XYZ', rg);
size(datRange)
siRange


%%
clc
rawDatRange = imfrmtDataSubset(rawDat, rawFrmt, rawRg);

imfrmtRangeSize(rawSi, rawFrmt, rawRg)
size(rawDatRange)


%%

datRange2 = imfrmtReshape(rawDatRange, rawFrmt, frmt, refrom, reto, rshi);
size(datRange2)

datRange - datRange2










%%












