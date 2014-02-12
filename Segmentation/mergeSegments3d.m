function segm = mergeSegments3d(label, param)
%
% segm = mergeSegments3d(label)
%
% description:
%    take 2d segmented 3d image and combine segments if overlap is large,
%    correct for undersegmentation
%
% input:
%    label   w x h x z matrix of label, each z plane contains labels of a 2d segmentation
%    param   parameter struct with
%
% output:
%    segm    w x h x z matrix of label for full 3d segmenation   
%
% reference:
%    LaTorre et al, 3D segmentations of neuronal nuclei from confocal microscope image stacks,
%    Frontiers in Neuroanatomy 2013
% 
% See also: 



% set of labels 

labs = unique(label(:));

%


end

