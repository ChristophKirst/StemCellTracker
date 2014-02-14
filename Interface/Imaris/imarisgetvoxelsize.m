function varargout = imarisgetvoxelsize(varargin)
%
% xyzsizes = imarisgetvoxelsize(imaris)
% [xsize, ysize, zsize] = imarisgetvoxelsize(imaris)
%
% description:
%   returns the voxel sizes 
%
% output:
%   xyzsizes    [xsize, ysize, zsize]
%   *size       Voxel size in direction *
%
% See also: imarisgetsize, imarisgetextend

[imaris, varargin, ~] = imarisvarargin(varargin);

dsize = imarisgetsize(imaris, varargin{:});
[smin, smax] = imarisgetextend(imaris, varargin{:});

vsize = (smin - smax) ./ dsize;

if nargout <= 1
    varargout{1} = vsize;
    
else
    varargout{1} = vsize(1);
    varargout{2} = vsize(2);
    varargout{3} = vsize(3);
    
end

end
