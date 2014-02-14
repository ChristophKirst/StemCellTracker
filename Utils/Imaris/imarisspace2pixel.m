function varargout = imarisspace2pixel(varargin)
%
% pixel = imarisspace2pixel(space)
% [pixelX, pixelY, pixelZ] = imarisspace2pixel(spaceX, spaceY, spaceZ)
% ... = imarisspace2pixel(object, ...)
% ... = imarisspace2pixel(imaris, ...)
% 
% description:
%     maps Imaris space coordinates to pixel coordinates using
%
% input:
%     space coordinates as [x;y;z] or x,y,z
%
% output:
%     pixel coordinates as [h,w,l] or h,w,l
%
% note:
%     Imaris stoares surfaces in space coordinates, here we convert these to pixel coordinates
%
% See also: imarispixel2space

[imaris, varargin, nargin] = imarisvarargin(varargin);

if isimaristype(imaris, varargin{1}, 'DataSet')
   dataset = varargin{1};
   varargin = varargin(2:end);
   nargin = length(varargin);
else
   dataset = imarisgetdataset(imaris);
end


if ~ismember(nargin, [1 3])
    error('imarisspace2pixel: expect [x;y;z] or x,y,z coordinates.');
end

  
% Check and get the inputs
if nargin == 1
    if size(varargin{1}, 2) ~= 3
        error('imarisspace2pixel: The input is expected to be an (N x 3) matrix.');
    end
    uPosX = varargin{1}(:, 1);
    uPosY = varargin{1}(:, 2);
    uPosZ = varargin{1}(:, 3);
    
elseif nargin == 3
    if all([size(varargin{1}, 1) > 1 size(varargin{1}, 2) > 1]) ...
       || all([size(varargin{2}, 1) > 1 size(varargin{2}, 2) > 1]) ...
       || all([size(varargin{3}, 1) > 1 size(varargin{3}, 2) > 1])
      error('imarisspace2pixel: The input is expected to be a vector.');
    end
    uPosX = varargin{1};
    uPosY = varargin{2};
    uPosZ = varargin{3};
end

% Get voxel sizes
voxelSizes = imarisgetvoxelsize(imaris, dataset);
if isempty(voxelSizes)
   error('imarisspace2pixel: cannot infer scaling.');
end

posX = (uPosX - dataset.GetExtendMinX()) / voxelSizes(1) + 0.5;
posY = (uPosY - dataset.GetExtendMinY()) / voxelSizes(2) + 0.5;
posZ = (uPosZ - dataset.GetExtendMinZ()) / voxelSizes(3) + 0.5;

if nargout == 0 || nargout == 1
    varargout{1} = [posX(:) posY(:) posZ(:)];
elseif nargout == 2
    varargout{1} = posX;
    varargout{2} = posY;
elseif nargout == 3
    varargout{1} = posX;
    varargout{2} = posY;
    varargout{3} = posZ;
else
    error('imarisspace2pixel: Bad number of output arguments.');
end

end
