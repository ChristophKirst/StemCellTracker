function varargout = imarispixel2space(varargin)
%
% pixel = imarispixel2space(space)
% [pixelX, pixelY, pixelZ] = imarispixel2space(spaceX, spaceY, spaceZ)
% ... = imarispixel2space(object, ...)
% ... = imarispixel2space(imaris, ...)
% 
% description:
%     maps Imaris pixel coordinates to space coordinates using
%
% input:
%     pixel coordinates as [h,w,l] or h,w,l
%
% output:
%     space coordinates as [x;y;z] or x,y,z
%
% note:
%     Imaris stores surfaces in space coordinates, here we convert these from pixel coordinates
%
% See also: imarisspace2pixel

[imaris, varargin, nargin] = imarisvarargin(varargin);

if isimaristype(imaris, varargin{1}, 'DataSet')
   dataset = varargin{1};
   varargin = varargin(2:end);
   nargin = length(varargin);
else
   dataset = imarisgetdataset(imaris);
end


if ~ismember(nargin, [1 3])
    error('imarispixel2space: expect [h;w;l] or h,w,l coordinates.');
end

  
% Check and get the inputs
if nargin == 1
    if size(varargin{1}, 2) ~= 3
        error('imarispixel2space: The input is expected to be an (N x 3) matrix.');
    end
    uPosX = varargin{1}(:, 1);
    uPosY = varargin{1}(:, 2);
    uPosZ = varargin{1}(:, 3);
    
elseif nargin == 3
    if all([size(varargin{1}, 1) > 1 size(varargin{1}, 2) > 1]) ...
       || all([size(varargin{2}, 1) > 1 size(varargin{2}, 2) > 1]) ...
       || all([size(varargin{3}, 1) > 1 size(varargin{3}, 2) > 1])
      error('imarispixel2space: The input is expected to be a vector.');
    end
    uPosX = varargin{1};
    uPosY = varargin{2};
    uPosZ = varargin{3};
end

% Get voxel sizes
voxelSizes = imarisgetvoxelsize(imaris, dataset);
if isempty(voxelSizes)
   error('imarispixel2space: cannot infer scaling.');
end

posX = (uPosX + 0.5) *  voxelSizes(1) + dataset.GetExtendMinX();
posY = (uPosY + 0.5) *  voxelSizes(2) + dataset.GetExtendMinY();
posZ = (uPosZ + 0.5) *  voxelSizes(3) + dataset.GetExtendMinZ();

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
    error('imarispixel2space: Bad number of output arguments.');
end

end
