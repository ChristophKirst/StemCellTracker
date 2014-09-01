function [shift, quality] = align2Images(img1, img2, varargin)
%
% [shift, quality] = align2Images(img1, img2, param)
%
% description: 
%    globally align two images using method specified in param
%    the function acts as a wrapper for the align2ImagesByXXX routines
% 
% input:
%    img1,img2  images to be aligned
%    param      (optional) parameter struct with entries:
%               .alignment       'Optimization', 'RMS', 'Correlation', 'Hugin' ('RMS")
%               parameters for methods align2ImagesByXXX
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
      error('align2Images: alignment %s not in %s', align, var2char(alignnames));
end

fun = str2func(['align2ImagesBy' algn]);

[shift, quality] = fun(img1, img2, param);

end