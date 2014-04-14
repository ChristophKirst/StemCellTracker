function imgroi = findRegionOfInfluence(imglab, param)
% 
% imgroi = findRegionOfInfluence(imglab, param)
%
% description:
%    for a labeld image imglab returns the regions of influence 
%    not exceeding a maximal distance specified in param
%
% input:
%    imglab    labeled image
%    param     (optional) parameter struct with entries
%              distance.max  maximal distance for influence zone (inf)
%
% output: 
%   imgroi     labeld imaged including the regoin of influecnes of the labels in imglab
%
% note:
%   labeled regoins should not touch 

if nargin < 2
   param = [];
end

distmax = getParameter(param, 'distance.max', inf);

dist = bwdist(imglab > 0);
imgroi = watershed(dist);

if distmax < inf
   mask = dist < distmax;
   imgroi = immask(imgroi, mask);
end

% make label identical to imglab:

pix1 = regionprops(imglab, 'PixelIdxList');
pix1 = cellfun(@first, {pix1.PixelIdxList});

pix2 = regionprops(imgroi, 'PixelIdxList');
pix2 = {pix2.PixelIdxList};

for i = 1:length(pix1)
   p1 = pix1(i);
   
   ll = imglab(p1);
   lr = imgroi(p1);

   if ll ~= lr
      imgroi(pix2{lr}) = ll;
   end
end

end
