classdef ROIRectangle < ROI
   %
   % ROIRectangle shape represents a 2D or 3D rectangle
   %
   % note: for convenice we 
   
   properties
      p1 = [0; 0];  % coordinate of lower left (bottom) (as row vector to allow for [obj.p1] concateations etc
      p2 = [0; 0];  % coordinate of upper right (top), p2 >= p1 
   end
   
   methods
      function obj = ROIRectangle(varargin)  % basic constructor
      %
      % Rectangle()
      % Rectangle(corner1, p2)
      %
      
         if nargin == 0
            obj.p1 = [0;0];
            obj.p2 = [0;0];
         elseif nargin == 1
            obj.fromArray(varargin{1});
         elseif nargin == 2
            obj.p1 = varargin{1};
            obj.p2 = varargin{2};
         else
            error('%s: invalid argument number %g for constructor', class(obj), nargin);
         end

         obj.initialize();
      end
      
      function initialize(obj)
         obj.p1 = obj.p1(:);  % make sure we have columns
         obj.p2 = obj.p2(:);  % make sure we have columns

         for i = 1:obj.dim    % order lower left and upper right !
            if obj.p1(i) > obj.p2(i)
               warning('%s: rectangle corners in dim %g were not ordered ', class(obj), i);
               c1 = obj.p1(i);
               obj.p1 = obj.p2; 
               obj.p2 = c1;
            end
         end     
      end
      
      function d = dim(obj)
         d = length(obj(1).p1);
      end

      function b = boundingbox(obj, varargin)
         b = obj.copy();  % rhe rectangle is its bbox
      end
      
      function v = volume(obj, varargin)
         v = prod([obj.p2]-[obj.p1], 1);
      end
      
      function a = toArray(obj)
         a = [obj.p1, obj.p2];
      end
      
      function obj = fromArray(obj, a)
         obj.p1 = a(:,1);
         obj.p2 = a(:,2);
      end

      function m = mask(obj, si)
         m = zeros(si);
         
         for d = obj.dim:-1:1
            coords{d} = obj.p1(d):obj.p2(d);
         end
         m(coords{:}) = 1;
      end
      
      function n = npixel(obj)
         n = prod([obj.p2]-[obj.p1]+1, 1);
      end
      
      function o = overlap(obj, roi)
         switch class(roi)
            case 'ROIRectangle'
               pp1 = max(a.p1, roi.p1);
               pp2 = min(a.p2, roi.p2);
               if any(pp1 > pp2)
                  o = ROIRectangle();
                  return
               end
               o = ROIRectangle(pp1, pp2);
            otherwise
               error('%s: overlap with % s not implemented!', class(obj), class(roi));
         end
      end
      
   end
   
end