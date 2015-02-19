function imwriteTIFFStack(img, filename, varargin)
%
% imwritestack(stack, filename, varargin)
%
% description:
%    writes image stack to file
%
% input:
%    stack     stack of images, either XYZ, XYCZ, XYZC (Z might be T)
%    filename  filename of saved images
%
% See also: imwrite, TIFF

isize = size(img);

if length(isize) == 3 % gray scale stack
   imwrite(img(:,:,1), filename, varargin{:}, 'WriteMode', 'overwrite');
   for i = 2:isize(3)
      imwrite(img(:,:,i), filename, varargin{:}, 'WriteMode', 'append');
   end
else % color image stack
   imwrite(img(:,:,:,1), filename, varargin{:}, 'WriteMode', 'overwrite');
   for i = 2:isize(4)
      imwrite(img(:,:,:,i), filename, varargin{:}, 'WriteMode', 'append');
   end
end
 
end


