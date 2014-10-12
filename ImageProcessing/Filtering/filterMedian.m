function out = filterMedian(image, ksize, padding)
%
% out = filterMedian(image, ksize, padding)
%
% description:
%     replaces each pixel by median of pixels within radius
%
% input:
%     image            image to filter
%     ksize            h x w (x l) size of the fitler
%     padding          (optional) specify the padding for boundary pixel
%
% output:
%     out              filtered image
%
% See also: filterMeanShift, filterBilateral

dim = ndims(image);

if dim < 2 || dim > 4
   error('filterMedian: expect 2d or 3d gray scale image!');
end

image = double(image);

if nargin < 2
   ksize = 3;
end
if length(ksize) < dim
   ksize = repmat(ksize(1), dim,1);
else
   ksize = ksize(1:dim);
end

if nargin < 3
   padding = 'replicate';
end


switch dim
   case 2
      % matlab is inconsistent as usual: replicate padding option not available -> use own median filter
      % out = medfilt2(image, ksize, padding);
          
      if nargin < 4
         chuncksize = 1;
      end      
      
      out = filterFunction(image, ksize, 'median', padding, chuncksize);
      

   case 3


      if nargin < 4
         chuncksize = 1;
      end

      out = filterFunction(image, ksize, 'median', padding, chuncksize);

end

end


