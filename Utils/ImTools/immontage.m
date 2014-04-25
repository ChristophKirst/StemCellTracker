function h = immontage(stack, normalize, varargin)
%
% image = immontage(stack, normalize, varargin)
%
% description:
%     montage of a image stack
%
% input:
%     stack      image stack
%     varargin   inputs to montage
%
% output:
%     h          handle
%
% See also: montage

if nargin < 2 || isempty(normalize)
   normalize = 1;
end

im = impqlpermute(stack, [], 'matlab');
if normalize
   im = double(im);
   im = im / max(im(:));
end
  
h = montage(im, varargin{:});

end
