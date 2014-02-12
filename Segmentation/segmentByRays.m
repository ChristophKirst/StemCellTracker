function label = segmentByRays(image, imagegrad, seeds, param)
%
% label = findRayRadius(image, imagegrad, seeds, param)
%
% description: 
%    starts rays from r0 and stops them if certain criteria are matched
%    returns a array labeling the segmetned identities of the different pixels 
%
% input:
%    image      instensity image to be segmented
%    imagegrad  gradinet of the intensity image
%    seeds      bw image of the size of image with white pixels indicating the sseds for the rays
%    param      struct with parameter entries as in findRayRadius and 
%               .join.overlap    % minimal overlap needed to join 
%               
%
% output:
%    label      sgements labeled with an index in each pixel, 0 = not segmented
%
% note: 
%    this is a wrapper for findRayRadius that returns the result as a labeled image
%    resolving conflicts due to overalps
%
% See also: findRayRadius

[~, ~, r0, polx, poly] = findRayRadius(image, imagegrad, seeds, param)
%






