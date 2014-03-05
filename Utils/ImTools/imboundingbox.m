function [minpos, maxpos] = imboundingbox(image)
%
% [minpos, maxpos] = imboundingbox(image)
%
% description:
%    finds the boundin box in pixel coordinates of the bw image
%
% input:
%    image  bw image
%
% output:
%    minpos  lower corner of bbox
%    maxpos  upper corner of bbox

indx =  find(image > 0);

ii = imind2sub(size(image), indx);

minpos = min(ii,[],1);
maxpos = max(ii,[],1);

end