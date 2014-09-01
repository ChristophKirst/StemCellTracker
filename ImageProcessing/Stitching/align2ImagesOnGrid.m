function [shift, quality] = align2ImagesOnGrid(imgs, varargin)
%
% [shift, quality] = align2ImagesOnGrid(imgs, varargin)
%
% desciption:
%    infers shifts of images in cell array imgs to make them align
%    the function is a wrapper that distributes the tasks to the the routines align2ImagesOnGridXXX
% 
% input:
%    imgs    images or image reading commands as cell array, the layout of the cell array usually should reflect 
%            the neighbouring relations between the images
%    param   (optional) parameter struct with entries:
%            .alignment       'Optimization', 'RMS', 'Correlation', 'Hugin'
%                             + parameters for methods align2ImagesOnGridByXXX
%                             ('RMS')
%
% output: 
%    shift   absolute shift of images assuming [lower,left(,bottom)] of image in cell array starts a [0,0(,0)] 
%            in pixel units and pixel coordinates
%    quality (optional) quality measure of alignment
%
% See also: alignImages


alignnames = {'Optimization', 'Correlation', 'RMS', 'Hugin'};
param = parseParameter(varargin{:});

algn = getParameter(param, 'alignment', 'RMS');
if ~any(strcmp(alignnames, algn))
      error('alignImagesOnGrid: alignment %s not in %s', align, var2char(alignnames));
end

fun = str2func(['align2ImagesOnGridBy' algn]);

[shift, quality] = fun(imgs, param);

end