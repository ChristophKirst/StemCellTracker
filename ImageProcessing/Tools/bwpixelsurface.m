function surface = bwpixelsurface(bw, method)
%
% surface = bwpixelsurface(bw)
%
% description: 
%    returns pixels on the surface of the bw objects
%    these are pixels with city-block distance 1 from the background = 0 
%
% input:
%    bw       the bw image
%    method   (optional) specify to include image borders or not ('border')
%             'border' or 'noborder'
%
% output:
%    surface  surface pixels of the bw region
%
% See also: impixelsurface, imsurface

if nargin < 2
   method = 'border';
end

switch method
   case 'border'
      %padd with background
      
      %ndims(bw)
      %size(bw)
      %ndims(bw)
      
      pad = padarray(bw, ones(ndims(bw), 1), 0);
      surface = bwdist(~pad, 'cityblock') == 1;
      surface = unpadarray(surface, ones(ndims(bw), 1));
      
   otherwise
      surface = bwdist(~bw, 'cityblock') == 1;
      %surface = bwdist(~bw, 'chessboard') == 1;
end

% alternative
%surface = imfilter(double(bw), ker); % imfilter uses 0 for padding by default
%surface(surface <= conn -2) = 0;
%surface(surface > 0) = 1;

end