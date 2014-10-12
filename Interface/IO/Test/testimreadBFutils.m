%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test various data access routines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

initialize
bfinitialize


%% Reader

ir = imreadBFReader('./Test/Images/hESCells_Colony.zvi')

%% Channel Names

cns = imreadBFChannelNames(ir, 'S', 1)
var2char(cns)

%% Positions

pp = imreadBFPositions(ir, 'S', [1,2,3], 'C', [1,2]);  
var2char(pp)

%% TimeStamps

ts = imreadBFTimestamps(ir);
var2char(ts)

%% Exposure Times

ts = imreadBFExposureTimes(ir);
var2char(ts)

%% Voxel Sizes (um)

vs = imreadBFVoxelSize(ir);
var2char(vs)


%% Wave Lengths (nm)

[em, ex] = imreadBFWavelengths(ir);
var2char(em)
var2char(ex)

%% Tags

tr = imreadBFTagRange(ir)

%%

tr = imreadBFTagRange(ir, 'y', 1:3, 'C', 1, 'X', [1,2,3])


