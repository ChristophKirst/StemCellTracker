function [ovl1, ovl2] = overlap2AlignedImages(img1, img2, sh)
%
% [ovl1, ovl2] = overlap2AlignedImages(img1, img2, sh)
%
% description:
%    returns the sub-images of img1 and img2 that overlap according to shift
%
% intput:
%    img1, img2   images
%    sh           shift between images in pql 
%
% output:
%    ovl1, ovl2   the overlapping regoing extracted form img1 and img2


s1 = size(img1); s2 = size(img2);
dim = ndims(sh);
sh0 = zeros(1,dim);
ov = findOverlap([sh0 + 1; s1], [sh + 1; sh + s2]);

if isempty(ov)
   ovl1 = [];
   ovl2 = [];
else
   ovl1 = imextract(img1, ov(1,:), ov(2,:));  
   ovl2 = imextract(img2, ov(1,:)-sh, ov(2,:)-sh);  
end

end

% helper   
function ov = findOverlap(a,b)
   a1 = a(1,:); a2 = a(2,:);
   b1 = b(1,:); b2 = b(2,:);
   
   ov = zeros(2,length(a1));
   ov(1,:) = max(a1,b1);
   ov(2,:) = min(a2,b2);
   
   if any(ov(2,:)-ov(1,:) < 0)
      ov = [];
   end
end