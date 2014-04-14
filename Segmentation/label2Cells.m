function [objects, varargout] = label2Cells(imglab, img, stats, param)
%
% objects = label2Cells(imglab, img, stats, param)
% [objects, stats] = label2Objects(...)
%
% description:
%    converts a labeld image and associated data to an array of Object classes
%    used for tracking and furtehr analysis
%
% input:
%    imglab  labeled image (2D/3D)
%    img     original grayscale image to measure basic properties from (2D/3D)
%    stats   (optional) previously calcualted statistics
%    param   (optional) parameter struct with entries:
%            .time      time for objects (0)
%            .rescale   rescale coordinates r by this factor ([1, 1(, 1)])
%            .method    how to calcualte the intensity in Object tracking, a string of any function, 'none' = dont calcualte ('median')
%            ,
%  
% output:
%    objects array of Objects each representing one of the labels
%    stats   (optional) updated statistics
%
% See also: Object, Cell, label2Cells

if nargin < 3
   stats = struc();
end

if nargin < 4
   param = [];
end

time = getParameter(param, 'time', 0);
if isempty(time)
   time = 0;
end

dim = ndims(imglab);

rescale = getParameter(param, 'rescale', 1);
if isempty(rescale)
   rescale = 1;
end
rescale = rescale(:)';
rescale = padright(rescale, dim, 'circular');

meth = getParameter(param, 'method', 'median');

stats = imstatistics(imglab, stats, {'Volume', 'Centroid', 'PixelIdxList'}, img);

if nargout > 1
   varargout{1} = stats;
end

nobj = length(prop);

switch meth
   case 'none'
      for i = nobj:-1:1
         objects(i) = Object(imglab(stats(i).PixelIdxList(1)), time, rescale' .* stats(i).Centroid', stats(i).Volume, 0);
      end
      
   otherwise
      fun = str2func(meth);
      
      for i = nobj:-1:1
         objects(i) = Object(imglab(stats(i).PixelIdxList(1)), time, rescale' .* stats(i).Centroid', stats(i).Volume, fun(image(stats(i).PixelIdxList)));
      end
end

end


