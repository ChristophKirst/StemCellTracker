function img = stitchImages(imgs, shifts, param)
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
% See also: alignImages, stitchImagesByMean, stitchImagesByMax, stitchImagesByMin, stitchImagesByOverwrite

if nargin < 3
   param = [];
end

meth = getParameter(param, 'method', 'Max');
if ~any(strcmp({'Mean', 'Max', 'Min', 'Overwrite'}, meth))
   error('stitchImages: method %s not Mean, Min, Max or Overwrite', meth);
end

fun = str2func(['stitchImagesBy' meth]);

img = fun(imgs, shifts);

end







