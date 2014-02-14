function [vertices, faces, normals] = imarisgetsurface(varargin)
%
% [vertices, faces, normals] = imarisgetsurface(timepoint)
% [vertices, faces, normals] = imarisgetsurface(objectname, ...)
% [vertices, faces, normals] = imarisgetsurface(object, ...)
% [vertices, faces, normals] = imarisgetsurface(imaris, ...)
%
% description:
%    get surface data from Imaris.
%
% input:
%    timepoint                 (optional) timepoint to put data (0)
%
%    objectename    name of Imaris surface
%    object         Imarise ISurface surface
%    imaris         Imaris application instance
%
% output:
%    vertices, faces, normals  surface data
%
% note: 
%    imaris returns suface in space coordinates, here reverted to double precision pixel coordinates
%    and convert faces indices to matlab format (+1)
%
% See also: imarisget

[imaris, varargin, nargin] = imarisvarargin(varargin);

if nargin < 1 || nargin > 2
   error('imarissetvolume: expect 1-2 input arguments');
end

if ischar(varargin{1})
   surfacename = varargin{1};
   surface = imarisgetobject(imaris, surfacename, 'Surface');
   if isempty(surface)
      error('imarisgetsurface: cannot find surface %s', surfacename);
   end 
   varargin = varargin(2:end);
   nargin = length(varargin);
elseif isimaristype(imaris, varargin{1}, 'Surface')
   surface = varargin{1};
   varargin = varargin(2:end);
   nargin = length(varargin);
else
   surface = imarisgetcurrentobject(imaris, 'Surface');
end


if nargin < 1
   timepoint = 0;
else
   timepoint = varargin{1};
end

nsurfaces = surface.GetNumberOfSurfaces;


ndata = 1;
vertices = {}; 
if nargout > 1
   faces = {}; 
end
if nargout > 2
   normals = {};
end


psize = imarisgetsize(imaris);
extend = imarisgetextend(imaris);
fac = psize./extend;

for i = 1:nsurfaces

      if surface.GetTimeIndex(i) == timepoint
         vertices{ndata} = surface.GetVertices(i); %#ok<AGROW>
         %vertices{ndata} = imarispace2pixel(imaris, vertices{ndata}); %#ok<AGROW>
         vertices{ndata} = imspace2pixel(psize, extend, vertices{ndata}); %#ok<AGROW>
      
         if nargout > 1
            faces{ndata} =  surface.GetTriangles(i) + 1; %#ok<AGROW>
         end
         if nargout > 2
            normals{ndata} = surface.GetNormals(i); %#ok<AGROW>
            %rescale 
            normals{ndata} = normals{ndata} * repmat(fac, size(normals{ndata},1),1); %#ok<AGROW>
         end
      end
      ndata = ndata + 1;
end


end
