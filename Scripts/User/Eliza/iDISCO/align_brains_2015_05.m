%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Align Brains from Imaris File %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initialize
ijinitialize

clc
close all

% make hdf5 8 work
setenv('HDF5_DISABLE_VERSION_CHECK', '2')

% make cmtk work
setenv('PATH', [getenv('PATH') ':' '/home/ckirst/Programs/cmtk-3.2.3/build/bin/'])


%%
basedir = '/home/ckirst/Science/Projects/BrainActivityMap/iDISCO_2015_05/';
fname1 = fullfile(basedir, 'brain 1/autofluo.h5');
fname2 = fullfile(basedir, 'brain 2/autofluo.h5');

if ~isfile(fname1)
   fprintf('file: %s does not exists!', fname1);
end


%%
fi = hdf5info(fname1)

%%

grps = fi.GroupHierarchy.Groups
{grps.Name}

%%

dname1 = '/t00000/s00/4/cells';
dname2 = '/t00000/s00/4/cells';

data1 = hdf5read(fname1, dname1);
size(data1)

data2 = hdf5read(fname2, dname2);
size(data2)

%%

%figure(5)
%implot(imclip(mat2gray(data), 0.05, 0.1))


%%
cllo = 0;
clhi = 0.25;
iju = ijplot3d(imclip(mat2gray(data1), cllo, clhi), 'Color', [1,0,0], 'Name', 'A');
iju = ijplot3d(imclip(mat2gray(data2), cllo, clhi), 'Universe', iju, 'Color', [0,0,1], 'Name', 'B')


%% Plot as Slices

figure(1); clf
implottiling(data1(:,:,1:6:end), 'tiling', [5,5])

figure(2); clf
implottiling(data2(:,:,1:6:end), 'tiling', [5,5])


%% Overlay

fig = figure(3); clf

data1RGB = imgray2color(data1(:,:,1:6:end), [1,0,0]);
data2RGB = imgray2color(data2(:,:,1:6:end), [0,1,1]);

implottiling(mat2gray(0.5* data1RGB + 0.5*data2RGB), 'tiling', [5,5])


%%

saveas(fig, fullfile(basedir, 'autofluo_1_2_overlay_raw.pdf'))


%% Write as Nrrd 

nrrdname1 =  fullfile(basedir, 'A.nrrd');
nrrdname2 =  fullfile(basedir, 'B.nrrd');

imwriteNRRD(nrrdname1, data1);
imwriteNRRD(nrrdname2, data2);

%% Run CMTK

wname = fullfile(basedir, 'warp.xform');

cmd = ['cmtk warp -o ', wname, ' --grid-spacing 100 --fast --jacobian-weight 1e-6 --refine 4 --verbose-level 9 --adaptive-fix-thresh 0.3 ', nrrdname1, ' ', nrrdname2];

cmd = ['cmtk warp -o ', wname, ' --grid-spacing 100 --fast --jacobian-weight 1e-6 --refine 4 --verbose-level 9 ', nrrdname1, ' ', nrrdname2];



cmd
%--rigidity-weight
%warp -o ffd5.xform --grid-spacing 40 --refine 3 --energy-weight 1e-1 \
%--initial affine.xform ref.nii flt.nii
%warp -o ffd5.xform --grid-spacing 40 --refine 3 --jacobian-weight 1e-5 \
%--initial affine.xform ref.nii flt.nii


%%

cmd = ['cmtk fview ', wname];
cmd


%% Reformat Volume

oname = fullfile(basedir, 'B_transform.nrrd');
cmd = ['cmtk reformatx -o ', oname, ' --floating ', nrrdname2, ' --verbose-level 9 ', nrrdname1, ' ', wname]

system(cmd)

%% Read Transformed Images an Plot

data1T = imreadNRRD(nrrdname1);
data2T = imreadNRRD(oname);

fig2 = figure(4); clf

data1RGB = imgray2color(data1T(:,:,1:5:end), [1,0,0]);
data2RGB = imgray2color(data2T(:,:,1:5:end), [0,1,1]);

implottiling(mat2gray(0.5* data1RGB + 0.5*data2RGB), 'tiling', [5,5])


%%

saveas(fig2, fullfile(basedir, 'A_B_overlay_aligned.pdf'))


