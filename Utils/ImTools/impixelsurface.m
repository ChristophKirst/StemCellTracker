function [pixsurf, varargout] = impixelsurface(imglab, stats)
%
% pixsurf = impixelsurface(imglab)
% [pixsurf, stats] = impixelsurface(imglab, stats)
%
% description: 
%    returns the pixel surface of the labeld image
%
% input:
%    label    the labeled image
%    stats    (optional) statistics of labeled image
%
% output:
%    pixsurf  labeled pixels on the pixsurf of the labeled objects
%    stats    (optional) updated statistics of labeld image
%
% See also: bwpixelsurface

if nargin < 2
   stats = [];
end

stats = imstatistics(imglab, stats, 'PixelIdxList', 'BoundingBox');

isize = size(imglab);
dim = length(isize);
if dim < 2 || dim > 3
   error('impixelsurface: expect 2d or 3d labeled image');
end

pixsurf = zeros(isize);

for l = 1:length(stats)
   if ~isempty(stats(l).PixelIdxList)
      bbox = stats(l).BoundingBox;
      
      %flat objects with trailling size 1 do not work with matlab !
      if bbox(dim)==bbox(2*dim)
         if bbox(2*dim) < isize(dim)
            bbox(2*dim) = bbox(dim) + 1;
         else
            bbox(dim) = bbox(dim) - 1;
            if bbox(dim) < 1
               bbox(dim) = 1;
            end
         end
      end
      
      bmin = bbox(1:dim);
      
      obj = imextract(imglab, bbox);
      obj = (obj == l);
      
      objsurf = bwpixelsurface(obj, 'border');
      idxsurf = find(objsurf);
      
      if dim == 2
         
         [sx,sy] = ind2sub(size(obj), idxsurf);
         sx = sx + bmin(1) - 1;
         sy = sy + bmin(2) - 1;
         pixsurf(sub2ind(isize, sx, sy)) = l;
         
      else
         
         [sx,sy,sz] = ind2sub(size(obj), idxsurf);
         sx = sx + bmin(1) - 1;
         sy = sy + bmin(2) - 1;
         sz = sz + bmin(3) - 1;
         pixsurf(sub2ind(isize, sx, sy, sz)) = l;
         
      end
   end
end

if nargout > 1
    varargout{1} = stats;
end

end