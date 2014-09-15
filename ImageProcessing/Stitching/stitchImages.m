function img = stitchImages(imgs, shifts, varargin)
%
% img = stitchImages(imgs, shifts, param)
%
% description:
%     stitches images imgs using shifts shifts and method specified in param
%
% input:
%     imgs    images to be stiched
%     shifts  shifts of images w.r.t. imgs{1}
%     param   (optional) parameter structu with entries:
%             .method 'Mean', 'Max', 'Min', 'Overwrite' ('Max')
%
% output:
%    img      stiched image
%
% See also: alignImages, stitchImagesByMean, stitchImagesByMax, stitchImagesByMin, stitchImagesByOverwrite, stitchImagesByHugin

param = parseParameter(varargin{:});

meth = getParameter(param, 'method', 'Max');
methnames = {'Mean', 'Max', 'Min', 'Overwrite', 'Hugin'};

if ~any(strcmp(methnames, meth))
   error('stitchImages: method %s not %s', meth, var2char(methnames));
end

fun = str2func(['stitchImagesBy' meth]);

img = fun(imgs, shifts, param);

end







