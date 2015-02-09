function [imgout, stats] = imlabelapplybw(imglab, fun, stats)
%
% labela = imlabelapplybw(imglab, fun)
%
% description:
%      applies fun that transfroms a bw image into a bw image onto each label separately
%
% input:
%      imglab    labeled image
%      fun       function operating on bw image representing a single region
% 
% output:
%      imglab    results of applying fun
%
% See also: imlabel

% lab = imlabel(label);
% labela = label;
% 
% for l = lab
%    imgl = fun(label == l);
%    labela(imgl > 0) = l;
% end

isize = size(imglab);
dim = length(isize);

if dim < 2 || dim > 3
   error('imlabelapplybw: expect 2d or 3d labeled image');
end

if ~exist('stats', 'var')
   stats = imstatistics(imglab, {'PixelIdxList', 'BoundingBox'});
else
   stats = imstatistics(imglab, stats, {'PixelIdxList', 'BoundingBox'});
end

n = length(stats);

imgout = zeros(isize);

for l = 1:n
   if mod(l, 500) == 0
      fprintf('imlabelapplybw: %g / %g\n', l, n);
   end
   
   idxpix = stats(l).PixelIdxList;
   
   if ~isempty(idxpix)
      ll = imglab(idxpix(1));
      
      bbox = stats(l).BoundingBox;
      
      %flat objects with trailing size 1 do not work with matlab !
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

      obj = imextract(imglab, bbox);
      obj = fun(obj == ll);

      idx = find(obj);
      
      if ~isempty(idx)
         sub = imind2sub(size(obj), idx);
         sub = sub + repmat(bbox(1:dim)'-1, size(sub,1), 1);
         ind = imsub2ind(isize, sub);
         imgout(ind) = ll;
      end
   end
end

end