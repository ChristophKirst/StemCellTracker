%%%% Test read of tracked supervoxels
% 

%% read the config file to set the various file names and path.

proj_dir = '/data/Science/Projects/StemCells/Experiment/Mouse/Nanog/TGMM/';
configFile = [proj_dir, 'TGMM_configFile.txt'];

% stores config file informtion internally, eqully well use
% configFile2FileNames(configFile, [], []);
configFile2FileNames(configFile);


%% read one 3D tiff

img_file = configFile2FileNames([], [], 5);
img = readTGMMStack(img_file, 1);
figure
implottiling(img(:,:,3:9), 'tiling', [3,2])

%% overlay supervox on image file, to check segmentation. The print below shows that contrary to TGMM docs that 
% some of the supervox from readListSupervoxelsFromBinaryFile() are of low quality and not
% in cells. This block allows one to look at single time data before computing
% trajectories

frame = 10;
UID = []; %'926';  % the last digits of an output/GMEM* files
img_file = configFile2FileNames([], UID, frame);
img = readTGMMStack(img_file, 1);
[cells, cell_label, sv_label, all_sv] = readOneTime(UID, frame, 1);
figure, hist([cells.Volume], 20);
title('histogram of cell volume (sum-z pixels)');
% uncomment for various statistics.
% stats = quickLabelImgStats(cell_label);
% printLabelImgStats(cell_label, stats);

% figure
% implottiling(sv_label(:,:,6:2:12), 'tiling', [2,2])
unique_sv = length(unique(sv_label(:))) - 1;  % -1 for background label=0
last_sv = max(all_sv(:));
fprintf('max sv label= %d, #sv-in-cells= %d\n', last_sv, unique_sv);

%% can play with various slices here. NB sv_label are sv in cells, all_sv is everything
% its very informative to offset by one frame the img data and the segmentation to see
% how much things move.
slice = 21;
[~, fname, ext] = fileparts(img_file);
imgrgb = overlayImgSuperVoxelBndry(img, all_sv, cell_label, slice);
figure, imshow(imgrgb)
title(sprintf('slice= %d, from file= %s. labeledImg: cells + SV ',slice, fname));

%% read all the times and assemble trajectories, 
% The numbering of frames in the tif file names are the last two arguments in the parse
% function, and define the matrix indices eg 1,2,.. in subsequent data files. See field
% all_frames_struct.fileId for how to recover the file name from matrix index.
%   One might save as .mat, the output of this cell, it can take 10 min for 100 frames

% the trajectory class should wrap up the statements in this block, stopping with the
% verifyTrajectory. The relevant objects (ie nuclei) are defined here via frames_struct
% and their statistics are computed and attached to the struct.
%UID = '280';
[all_frames_struct, trajectory] = parseMixtureGaussians2FrameTrajectory(UID, 0, 15);
% add various stats to trajectory that are inherited from cells. Do it before filtering out bad
% trajectories, since minimal trajectories mark births and should be retained to aid visualization.
trajectory = addRegionprops2Trajectory(all_frames_struct, trajectory);
quickTrajectoryStats(trajectory);
% code to add stats to cells/nuclei as listed in frames, which can then be used to filter
% trajectories in verifyTrajectories. Can then do graphics based on
% trajectories. When add other data REREAD cell pixel list from binary files and do
% averages, do not save pixelIdxList. 
trajectory = verifyTrajectory(all_frames_struct, trajectory);
% select trajectories that pass filter, ie quality = 1
good_traj = trajectory( logical([trajectory.quality]) );
% optional save incase need restart due to ijViewCZT crashing Java.
%save_fn = [proj_dir filesep 'traj' UID '.mat'];
%save(save_fn, 'all_frames_struct', 'trajectory');

%% overlay other color channels on the nuclear data and compute statistics
% Use frames and trajectories from last block

% in_dir = '/Users/siggia/test';
% pattern = 'ColonyW02T_T<4>_C2.tif';
% all_frames_struct = addChannelStats(UID, all_frames_struct, in_dir, pattern, 1:10);
%% visualize a C(olor)ZT array with img overlayed with cell boundaries labeled with trajectory.Id
% be aware of the final image sizes. Must run previous block first..Can also restrict
% attention here to just 'good' trajectories or anything else.
%   To diagnose why certain cells are not tracked, use 
% label_img = trajectory2LabelImg(UID, all_frames_struct, trajectory, file_frames, 'filter', 0);
% which does not omit anything, and you will discover how the tracking is going wrong.

%UID = '249'; %'294'; %
file_frames = 1:5; % NB frame numbers in file names, can begin with 0 or 1
% create a labeled image based on trajectory.id x,y,z,times. The options on following
% routine restrict to trajectories that span the frame interval. Use filter,0 option to
% see all trajectories (and almost all cells) to understand why the tracking is off
label_img = trajectory2LabelImg(UID, all_frames_struct, trajectory, file_frames, 'filter', 0, 'family', 0);
% add the nuclear label density and color code cell boundaries.
imgCZT = trajectory2IJ(UID, label_img, file_frames);
[~, ~, nz, nt] = size(label_img);
% view the whole thing, if not too large  if not eg imgCZT(100:600,100:600,:)

imgCZT = reshape(imgCZT, size(imgCZT,1), size(imgCZT,2), 3, nz, nt);
size(imgCZT)

imgZCT = imfrmtReformat(imgCZT, 'XYCZT', 'XYZCT');
ijplot5d(imgZCT)

%ijViewCZT(imgCZT(:, :,:), nz, nt, sprintf('frames= %d..%d', file_frames(1), file_frames(end) ))

%% cutout region around selected trajectory in prev big images and plot. Selected traj's
% highlighted with white boundary. Redo the labeled image to include all trajectories so
% can see all the births, even resulting in bad trajectories ie too short.
file_frames = 1:10;
label_img = trajectory2LabelImg(UID, all_frames_struct, trajectory, file_frames, 'filter', 0);
imgCZT = trajectory2IJ(UID, label_img, file_frames);

[img_czt, nz, ntimes] = highlightTrajectory(label_img, imgCZT, good_traj(1) );


%ijViewCZT(img_czt, nz, ntimes, 'test of cutout');

img_zct = imfrmtReformat(img_czt, 'XYCZT', 'XYZCT');
ijplot5d(img_zct)

% check out the testExamineTrajectory script
%%%%%%%%%%%%%%%%%%%% rest not of general interest %%%%%%%%%%%%%%%%%
%% lower level reads of data. The cells struct array here, is lacking some stats supplied by readOneTime.

[~, data_file] = configFile2FileNames([],[],5);
[svList, sizeIm] = readListSupervoxelsFromBinaryFile( [data_file '.svb'] );
cells = xmlMixtureGaussians2struct( [data_file '.xml'] ); 
all_pix = sum( cellfun(@length, svList));
ncells = length(cells);
fprintf('ncells= %d, #super-vox= %d, total vox= %d pct-vol= %3.2f, sizeImg= %d %d %d\n',...
    ncells, length(svList), all_pix, 100*all_pix/prod(sizeIm), sizeIm);


%% convert the svList for one super-vox into pixel coordinates of the super-vox.

[~, data_file] = configFile2FileNames([],[],5);
[svList, sizeIm] = readListSupervoxelsFromBinaryFile( [data_file '.svb'] );
ptr_sv = 1;
[i,j,k] = ind2sub(sizeIm, svList{ptr_sv});
fprintf('len supervoxel= %d mean,std i,j,k %3.1f +- %3.1f, %3.1f +- %3.1f, %3.1f +- %3.1f\n',...
    length(i), mean(i),std(i), mean(j),std(j), mean(k),std(k) );

%% routine from TGMM, that processes all times, and computes phylogenies of cells from
% parseMixtureGaussiansXml2trackingMatrixCATMAIDformat.m
% my file parseMixtureGaussians2FrameTrajectory has the same functionality
base = '/Users/siggia/Desktop/TGMM/H2BCitrineTest/output/GMEMtracking3D_1417552310/XML_finalResult_lht/GMEMfinalResult_frame';
[trackingMatrix, svIdxCell] = parseMixtureGaussiansXml2trackingMatrixCATMAIDformat(base, 1, 20);

% svIdxCell has cell number for each of the supervoxels that are tracked.
allIdx = [svIdxCell{:}];

fprintf('\nsize svIdxCell= %d %d, len allIdx= %d, min,max super-voxIdx= %d %d\n',...
    size(svIdxCell), length(allIdx), min(allIdx), max(allIdx));