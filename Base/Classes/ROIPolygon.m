classdef ROIPolygon < ROI
   %
   % ROIPolygon shape represents a 2D polygon
   %
   % note: the coordinates of the polygon are stored as rows in .p array
   
   properties
      p = [];  % coordinates as rows
   end
   
   methods
      function obj = ROIPolygon(varargin)  % basic constructor
         %
         % Rectangle()         empty rectangle
         % Rectangle(p1, p2)   rectangle with opposite corners p1 and p2
         %
         
         if nargin == 0
            obj.p = [];
         elseif nargin == 1
            obj.fromPixelArray(varargin{1});
         else
            error('%s: invalid argument number %g for constructor', class(obj), nargin);
         end
      end

      function a = toArray(obj)
         a = obj.p;
      end
      
      function obj = fromArray(obj, a)
         obj.p = a;    
      end
      
      function obj = fromPixelArray(obj, a)  
         obj.p = a;
      end
      
      function a = toPixelArray(obj)
         a = round(obj.p);
      end

      
      function d = dim(obj)
         d = size(obj(1).p, 1);
      end

      function b = boundingBox(obj, varargin)
         b = ROIRectangle(min(obj.p,[],2), max(obj.p,[],2));
      end
      
      function v = volume(obj, varargin)
         v = polyarea(obj.p(1,:), obj.p(2,:));
      end
 
      function p = pixelIdxList(obj, si)
         p = find(obj.mask(si));
      end
      
      function m = mask(obj, si)
         x = obj.p(1,:);
         y = obj.p(2,:);
         m = poly2mask(y,x, si(1), si(2));
         m = imdilate(m, strel('square', 3)); % poly2 mask is very restrictive
         n = size(obj.p,2);
         %add lines
         %for i = 1:n
         %   m = impixelline(m, obj.p(:,i), obj.p(:,mod(i,n)+1), 1);
         %end
         x = round(x); y = round(y);
         ids = and(and(x > 0, x <= si(1)), and(y > 0, y <= si(2)));
         x = x(ids); y = y(ids);
         ids = sub2ind(si, x,y); % add corners
         m(ids) = 1;
      end
      
      function n = npixel(obj, si)
         n = total(obj.mask(si));
      end
      
      function o = overlap(obj, roi) %#ok<STOUT>
         error('%s: overlap with % s not implemented!', class(obj), class(roi));
      end
      
      
      % exract roi from an array / image
      function [d, sh] = extractData(obj, d)
         %
         % [d, sh] = extractdata(obj, d)
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
         
         
         bb = obj.boundingBox;
         [d, sh] = bb.extractData(d);
         shr = repmat(sh(:), 1, size(obj.p,2));
         obj.p = obj.p - shr;
         m = obj.mask(size(d));
         obj.p = obj.p + shr;
         d = immask(d, m);
      end

      function shift(obj, sh)
         obj.p = obj.p + repmat(sh(:), 1, size(obj.p,2));
      end

      
      function obj= dilate(obj, bwidth)
         obj.p = dilatePolygon(obj.p, bwidth);
      end
      
      
      function plot(obj, varargin)
         n = length(obj);
         cc = colorcube(n);
         for i = 1:n
            pp = obj(i).p;
            plot(pp(1,:)', pp(2,:)', varargin{:}, 'Color', cc(i,:))
         end
      end
      
      function obj = rescale(obj, scale)
         for i = 1:length(obj)
            obj(i).p = obj(i).p .* scale;
         end
      end
  
   end
   
end