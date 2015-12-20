%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test Writing Points to Imaris File %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initialize
ijinitialize

clc
close all

% make hdf5 8 work
setenv('HDF5_DISABLE_VERSION_CHECK', '2')


%%

fn = '/home/ckirst/Desktop/test for spots without spot.ims'
isfile(fn)

%%
fi = hdf5info(fn)

%%

fi.GroupHierarchy

%% Create Points

pts = 1:10;
npts = length(pts);
t = 0 .* pts;
pts = [1900 + 50 * sin(2 * pi * pts / npts); 2000 + 50 * cos(2 * pi * pts / npts); 1700 + 0 .* pts; 2.5 + 0 .* pts]'


%%
pn = '/Scene/Content/Points0'

pnt = [pn, '/Time']
%hdf5write(fn, pnt, [npts, 1], 'DataType', 'int64');
hdf5write(fn, pnt, int64(t));



