function threshold = otsuThreshold(image)
%
% threshold = otsuThreshold( image )
%
% description:
%    finds Otsu threshold of the image
%
% input:
%    image       the image to calculate threshold for
%
% output:
%    threshold   calculated threshold
%
% See also: graylevel

threshold = graythresh(image);

end

