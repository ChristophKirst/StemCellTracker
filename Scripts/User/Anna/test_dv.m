
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Serial Detection of Colonies 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% image data is read from disk each time and not cached
% useful for large tilings that do not fit into memory

clear all
clear classes
close all
clc
initialize
bfinitialize

verbose = true;

initializeParallelProcessing(8)

%%

fn = '/home/ckirst/Science/Projects/StemCells/Data/TrackingPC/Anna/11082014/LY40_DAPI_OCT4_NANOG_SOX2_01_R3D.dv';
is = ImageSourceBF(fn);
is.setCellDataFormat('XY', 'SC')
is.addRange('C', 1)
is.printInfo

%% Stage Positions

pos = is.stagePositions('C', 1);
pos = pos(:);

vsiz = is.voxelSize('S', 1)
vsiz = vsiz{1}

figure(1); clf
ppos = [pos{:}];
plot(ppos(1,:)', ppos(2,:)')

%% Alignment

clc
pos2d = cellfunc(@(x) x(1:2), pos);

algn = Alignment;
algn.asource = is;
algn.anodes = 1:is.nCells;

algn.fromStagePositionsAndVoxelSize(pos2d, vsiz)

%%

figure(3)
ppos2d = [pos2d{:}];
plot(ppos2d(1,:)', ppos2d(2,:)', '*')

%%

posI = algn.imagePositions;
posI = cellfunc(@(x) x', posI);
figure(4)
pposI2d = [posI{:}];
plot(pposI2d(1,:)', pposI2d(2,:)', '*')

%%

figure(2);
algn.plotPreviewStiched

%%




