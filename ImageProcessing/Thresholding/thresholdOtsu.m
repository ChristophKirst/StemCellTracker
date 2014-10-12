function threshold = thresholdOtsu(image)
%
% threshold = thresholdOtsu( image )
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

