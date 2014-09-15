classdef Disk < Shape
   %
   % Disk shape represents a 2D or 3D disk / ellipsoid
   %
   
   properties
      radii = [0,0];
      center = [0,0];
      phi = 0;
   end
   
   methods
      function obj = Disk(varargin)  % basic constructor
      %
      % Rectangle()
      % Rectangle(center, radii)
      %
         if nargin > 0
            obj.center = varargin{2};
         end
         if nargin > 1
            obj.radii = varargin{1};
         end
         
         obj.initialize();
      end
      
      function initialize(obj)
         obj.radii = obj.radii(:);
         if  length(obj.radii) < obj.dim
            obj.radii = padright(obj.radii, obj,dim, obj.radii);
         end
         
         obj.center = obj.center(:);
      end
      
      function d = dim(obj)
         d = length(obj(1).center);
      end

   %% Todo:
      function b = boundingbox(obj, varargin)
         b = [, obj.corner2];
      end
      
      function v = volume(obj, varargin)
         v = prod(obj.corner2-obj.corner);
      end
      
   end
end
