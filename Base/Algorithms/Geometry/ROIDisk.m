classdef ROIDisk < Shape
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
            obj.center = varargin{2};
         end
         if nargin > 1
            obj.radius = varargin{1};
         end
         
         obj.initialize();
      end
      
      function initialize(obj)
         obj.radius = obj.radius(:); 
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
         for d = obj.dim:-1:1
            coords{d} = floor(obj.p1(d)+1):floor(obj.p2(d));
         end
         p = sub2ind(si, coords{:});
      end
      
      function m = mask(obj, si)
         m = zeros(si);
         m(obj.pixelIdxList(si)) = 1;
      end
      
      function n = npixel(obj)
         n = prod(floor([obj.p2])-floor([obj.p1])-1, 1);
      end
      
      function o = overlap(obj, roi)
         switch class(roi)
            case 'ROIRectangle'
               pp1 = max(a.p1, roi.p1);
               pp2 = min(a.p2, roi.p2);
               if any(pp1 > pp2)
                  o = ROIRectangle();
               else
                  o = ROIRectangle(pp1, pp2);
               end
            otherwise
               error('%s: overlap with % s not implemented!', class(obj), class(roi));
         end
      end
      
      
      
      
      function su = pixelsurfaceIdxList(obj, si)

            x = radius;
            y = 0;
            radiusError = 1-x;
            
            while (x >= y)
            
               sub2ind(x + x0, y + y0);
               DrawPixel(y + x0, x + y0);
               DrawPixel(-x + x0, y + y0);
               DrawPixel(-y + x0, x + y0);
               DrawPixel(-x + x0, -y + y0);
               DrawPixel(-y + x0, -x + y0);
               DrawPixel(x + x0, -y + y0);
               DrawPixel(y + x0, -x + y0);
               y++;
               if (radiusError<0)
               {
                  radiusError += 2 * y + 1;
                  }
               else
                  {
                     x--;
                     radiusError += 2 * (y - x + 1);
                     }
                  }
                  }
                  
         
      end
      
      
      
      % special functions for rectangle
      function pp = p(obj)
         pp = obj.p1(1);
      end
      
      function qq = q(obj)
         qq = obj.p1(2);
      end
      
      function ll = l(obj)
         ll = obj.p1(3);
      end
      
      function ww = w(obj)
         ww = obj.p2(1) - obj.p1(1);
      end
      
      function hh = h(obj)
         hh = obj.p2(2) - obj.p1(2);
      end
         
      function dd = d(obj)
         dd = obj.p2(3) - obj.p1(3);
      end
      
      
   end
      
      
  
end
