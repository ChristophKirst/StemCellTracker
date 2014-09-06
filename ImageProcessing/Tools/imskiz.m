function skiz = imskiz(image, label, mask)
%
% skiz = imskiz(image, seeds, mask)
%
% description: 
%     calculate the geodesic SKIZ (sceleton of influence zone) of the labels
%
% input:
%     image   input bw image
%     label   seed for the skizs
%     mask    (optional) mask for roi
%
% ouput: 
%     skiz    skeleton of influence zones as labeled image
%
% See also: bwdistgeodesic, watershed

% the skiz is given by the watershed of the geodesic distance transform

bwg = mat2gray(bwdistgeodesic(image, label > 0));

skiz = imwatershedlabel(bwg, label);
skiz = skiz .* cast(image, class(skiz));

if nargin > 2
   skiz = skiz .* cast(mask, class(skiz));
end

end

   


