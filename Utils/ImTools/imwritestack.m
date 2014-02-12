function imwritestack(stack, filename, varargin)
%
% imwritestack(stack, filename, varargin)
%
% description:
%    writes image stack to file
%
% input:
%    stack     stack of images
%    filename  filename of saved images
%    varargin  args to imwrite
%
% See also: imwrite

dim = size(img);

if length(dim) == 3 % gray scale stack
   imwrite(stack(:,:,1), filename, varargin{:}, 'WriteMode', 'overwrite');
   for i = 2:dim(3)
      imwrite(stack(:,:,i), filename, varargin{:}, 'WriteMode', 'append');
   end
else % color image stack

   imwrite(stack(:,:,:,1), filename, varargin{:}, 'WriteMode', 'overwrite');
   for i = 2:dim(4)
      imwrite(stack(:,:,:,i), filename, varargin{:}, 'WriteMode', 'append');
   end
end
 
end