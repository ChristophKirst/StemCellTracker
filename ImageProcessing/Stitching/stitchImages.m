function img = stitchImages(imgs, ipos, varargin)
%
% img = stitchImages(imgs,  isize, ipos, param)
%
% description:
%     stitches images imgs using shifts shifts and method specified in param
%
% input:
%     imgs    images to be stiched
%     ipos    image positions or shifts
%     param   (optional) parameter structu with entries:
%             .method 'Mean', 'Max', 'Min', 'Overwrite' ('Max')
%             .size   final image size ([] = automatic, positions are interpreted as shifts and converted to absolute shifts)
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

img = fun(imgs, ipos, param);

end







