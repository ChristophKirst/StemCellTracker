classdef ROIRectangle < ROI
   %
   % ROIRectangle shape represents a 2D rectangle or 3d cube
   %
   % note: the pixel of a rectangle from (0,0) -> (1,1)  is pixel (1,1) etc
   %       in particular a rectangle (0,0) -> (0,0) is empty and (0,0) -> isize spans an image fo size isize
   %       all coords arre assumed to be row vectors to allow for [obj.p1] concateations etc
   
   properties
      p1 = [0; 0];  % coordinate of lower left (bottom) 
      p2 = [0; 0];  % coordinate of upper right (top), p2 >= p1 
   end
   
   methods
      function obj = ROIRectangle(varargin)  % basic constructor
         %
         % Rectangle()         empty rectangle
         % Rectangle(p1, p2)   rectangle with opposite corners p1 and p2
         %
         
         if nargin == 0
            obj.p1 = [0;0];
            obj.p2 = [0;0];
         elseif nargin == 1
            obj.fromPixelArray(varargin{1});
         elseif nargin == 2
            obj.fromPixel(varargin{1}, varargin{2});
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

      function b = boundingBox(obj, varargin)
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
      
      function obj = fromPixel(obj, p1, p2)  
         obj.p1 = p1 - 1;
         obj.p2 = p2;
      end
      
      function obj = fromPixelArray(obj, a)  
         obj.fromPixel(a(:,1), a(:,2));
      end
      
      function a = toPixelArray(obj)
         a = round([obj.p1+1, obj.p2]);
      end

      function p = pixelIdxList(obj, si)
         p = find(obj.mask(si));
      end
      
      function m = mask(obj, si)
         m = zeros(si);
         rng = cell(1, obj.dim);
         l = round(max(obj.p1+1, 1));
         u = round(min(obj.p2, si(:)));
         for d = 1:obj.dim
            rng{d} = l(d):u(d);
         end
         m(rng{:}) = 1;
      end
      
      function n = nPixel(obj, si)
         if nargin < 2
            si = obj.p2;
         end
         n = prod(round(min(obj.p2, si(:)))-round(max(obj.p1, 1)), 1);
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
      
      
      % extract roi from an array / image
      function [d, sh] = extractData(obj, d)
         %
         % [d, sh] = extractData(obj, d)
         %
         % description:
         %     extracts data from boudnign box
         %
         % input:
         %     d    data
         %
         % output:
         %     d    extracted data
         %    sh    (optional) shift of lower left corner w.r.t to full image
         
         
         dim = ndims(d);
         if dim ~= obj.dim
            error('ROIRectangle: extractdata: data dimension mismatch');
         end
         rectlow = round(obj.p1 + 1);
         recthigh = round(obj.p2);
         si = size(d);
         sh = zeros(1,dim);
         for i = dim:-1:1
            sh(i)   = max(rectlow(i),1);
            rect{i} = sh(i):min(recthigh(i), si(i));
         end  
         d = d(rect{:});
         
      end


      function shift(obj, sh)
         obj.p1 = obj.p1 + sh(:);
         obj.p2 = obj.p2 + sh(:);
      end
      
      function obj = rescale(obj, scale)
         for i = 1:length(obj)
            obj(i).p1 = obj(i).p1 .* scale;
            obj(i).p2 = obj(i).p2 .* scale;
         end
      end
      
      
      function plot(obj, varargin)
         n = length(obj);
         cc = colorcube(n);
         for i = 1:n
            p11 = obj(i).p1; p22 = obj(i).p2;
            plot([p11(1), p11(1), p22(1), p22(1), p11(1)], [p11(2), p22(2), p22(2), p11(2),p11(2)] , varargin{:}, 'Color', cc(i,:))
         end
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