function imgseg = segmentOnTiling(image, segmenter, param)
%
% imgseg = segmentOnTiling(image, segmenter, param)
%
% description:
%    uses the segmentation routine segmenter to segment on a finer tiling 
%    for paralleization and processing speed up
%
% input:
%    image     3d grayscale or color image
%    param     (optional) parameter struct with entries
%              .color.scale    scale color data (true)
%              .color.alpha    alpha data (image intensity)
%              .renderslice    set the slices to render ('z')
%              .range.h        1x2 h-axis bounds. ([1 size(data, 1)])
%              .range.w        1x2 w-axis bounds. ([1 size(data, 2)])
%              .range.l        1x2 z-axis bounds. ([1 size(data, 3)])
%              .boxratios      1x3 ratio array to scale the axes individually ([1 1 1])
%              .parent         parent axes (gca)
%
% output:
%    model     refernece to the 3d model
%
% note:
%    use interp3 on input date to increase/decrease resolution of data
%
% See also alphamap, colormap, opengl, isosurface