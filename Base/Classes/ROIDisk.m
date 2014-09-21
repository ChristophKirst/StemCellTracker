classdef ROIDisk < ROI
   %
   % Disk shape represents a 2D disk or 3D sphere
   %
   
   properties
      radius = 0;
      center = [0,0];
   end
   
   methods
      function obj = ROIDisk(varargin)  % basic constructor
      %
      % Rectangle()
      % Rectangle(center, radii)
      %
         if nargin > 0
            obj.center = varargin{1};
         end
         if nargin > 1
            obj.radius = varargin{2};
         end
         
         obj.initialize();
      end
      
      function initialize(obj)
         obj.radius = obj.radius(1); 
         obj.center = obj.center(:);
      end
      
      function d = dim(obj)
         d = length(obj(1).center);
      end

      function v = volume(obj)
         if obj.dim == 2
            v = pi * [obj.radius].^2;
         else
            v = 4/3 * pi * [obj.radius].^3;
         end
      end

      function b = boundingbox(obj, varargin)
         b = Rectangle(obj.center - r, obj.center + r);
      end

      
      function p = pixelIdxList(obj, si)
         p = find(obj.mask(si));
      end
      
      function m = mask(obj, si)
         if nargin < 2
            si = obj.center + obj.radius;
         end
         m = imgrid(si);
         x = m{1}; y = m{2}; % assume 2d
         x0 = obj.center(1); y0 = obj.center(2);
         r = obj.radius;
         m = ((x- x0).^2 + (y-y0).^2 - r^2 <= eps);
      end
      
      function n = npixel(obj)
         n = total(obj.mask);
      end
      
      function o = overlap(obj, roi) %#ok<STOUT>
         error('%s: overlap with % s not implemented!', class(obj), class(roi));
      end

      function su = pixelsurfaceIdxList(obj, si)
         x0 = obj.center(1);
         y0 = obj.center(2);
         
         x = obj.radius;
         y = 0;
         radiusError = 1-x;
            
         su = zeros(2, ceil(2 * pi * obj.radius));
         i = 0;
         
         while (x >=y)
            su(:,i)   = [x + x0; y + y0];
            su(:,i+1) = [y + x0; x + y0];
            su(:,i+2) = [-x + x0; y + y0];
            su(:,i+3) = [-y + x0; x + y0];
            su(:,i+4) = [-x + x0; -y + y0];
            su(:,i+5) = [-y + x0; -x + y0];
            su(:,i+6) = [x + x0; -y + y0];
            su(:,i+7) = [y + x0; -x + y0];
            i = i + 8;
            y = y + 1;
            if (radiusError<0)
               radiusError = radiusError + 2 * y + 1;
            else
               x = x-1;
               radiusError = radiusError + 2 * (y - x + 1);
            end
         end
         
         su = su(:,1:i-8);
         su(su(1,:) <= 0 || su(1,:) > si(1), :) = [];
         su(su(2,:) <= 0 || su(2,:) > si(2), :) = [];
         su = imsub2ind(si, su);
      end
      

      function d = extractdata(obj, d)
         d = immask(d, obj.mask);
      end
      
      
      % special functions for disk
      function rr = r(obj)
         rr = obj.radius;
      end
      
      function pp = p(obj)
         pp = obj.center(1);
      end
      
      function qq = q(obj)
         qq = obj.center(2);
      end
      
      function ll = l(obj)
         ll = obj.center(3);
      end
      
      
   end
      
      
  
end
