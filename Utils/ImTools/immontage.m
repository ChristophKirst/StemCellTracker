function h = immontage(stack, varargin)
%
% image = immontage(stack, varargin)
%
% description:
%     montage on image stack
%
% input:
%     stack      h x w x z stack of h x w images
%     varargin   inputs to montage
%
% output:
%     h          handle
%
% See also: montage

if ndims(stack) ~= 3
   error('immontage: expect stack of grayscale images')
end

h = montage(imgray3d(stack), varargin{:});

end
