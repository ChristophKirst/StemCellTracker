function stats = overlapStatistics2AlignedImages(img1, img2, sh, varargin)
%
% stats = overlapStatisticsImagePair(img1, img2, sh)
%
% description:
%   this routine determines statistics in the overlap region of two aligned images
%
% input:
%    img1, img2   images
%    sh           shift between images
%
% output:
%    stats  struct with entries .from and .to
%           each having entries .var, .max, .min of potential overlap regoins
%

s1 = size(img1); s2 = size(img2);
dim = ndims(sh);
sh0 = zeros(1,dim);
ov = findOverlap([sh0 + 1; s1], [sh + 1; sh + s2]);

if isempty(ov)
   fr.var = 0;
   fr.max = Inf;
   fr.min = -Inf;
   
   stats.from = fr; 
   stats.to   = fr;
   return
end

% stats of img1 overlap region

img1 = imextract(img1, ov(1,:), ov(2,:));  
fr.var = var(double(img1(:)));
fr.max = double(max(img1(:)));
fr.min = double(min(img1(:)));


% stats of img2 overlap region

img2 = imextract(img2, ov(1,:)-sh, ov(2,:)-sh);  
to.var = var(double(img2(:)));
to.max = double(max(img2(:)));
to.min = double(min(img2(:)));

stats.from = fr;
stats.to   = to;

end


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


