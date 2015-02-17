classdef ROIPolygon < ROI
   %
   % ROIPolygon shape represents a 2D polygon
   %
   % Polygon is represented as a cell array of contours (assuming even odd order)
   % each contour is represented as 2xn matrix of coordinates (coorinates as rows)
   %
   % TODO: clean up and integrate fully to polygon package!
   
   properties
      contours = {};
   end
   
   methods
      function obj = ROIPolygon(varargin)  % basic constructor
         %
         % ROIPolygon()    empty polygon
         % ROIPolygon(p)   single contour polygon, p coordinates as rows
         % ROPPolygon(c)   c cell array of contours
         
         if nargin == 0
            obj.contours = {};
         elseif nargin == 1
            if isequal(class(varargin{1}), 'triangulation')
               obj.fromTriangulation(varargin{:});
            elseif iscell(varargin{1})
               obj.fromCell(varargin{:});
            else
               obj.fromPixelArray(varargin{1});
            end
         else
            error('%s: invalid argument number %g for constructor', class(obj), nargin);
         end
      end

      function a = toArray(obj)
         if length(obj.contours) == 1
            a = obj.contours{1};
         else
            error('ROIPolygon: toArray: polygon consists of multiple contours!');
         end
 
      end
      
      function obj = fromArray(obj, a)
         obj.contours = {a};    
      end
      
      function obj = fromCell(obj, c)
         obj.contours = c;
      end
      
      function obj = fromPolygon(obj, p)
         obj.contours = p;
      end
      
      function obj = fromPixelArray(obj, a)  
         obj.contours = {a};
      end
      
      function a = toPixelArray(obj)
         a = round(obj.toArray());
      end
      
      function p = toPolygon(obj)
         p = obj;
      end
      
      function t = toTriangulation(obj, varargin)
         t = polygonToTriangulation(obj.contours, varargin{:});
      end
      
      function obj = fromTriangulation(tri, varargin)
         obj.contours = triangulationToPolygon(tri, varargin{:});
      end

      
      function d = dim(~)
         %d = size(obj(1).contours{1}, 1);
         d = 2;
      end

      function b = boundingBox(obj, varargin)
         b = ROIRectangle(min(cellfun(@(x) min(x,[],2), obj.contours)), max(cellfun(@(x) max(x,[],2), obj.contours)));
      end
      
      function v = volume(obj, varargin)
         v = zeros(1,length(obj));
         for i = 1:length(obj)
            v(i) = polygonArea(obj.countours);
         end
      end
 
      function p = pixelIdxList(obj, si)
         p = find(obj.mask(si));
      end
      
      function m = mask(obj, si)
%          x = obj.p(1,:);
%          y = obj.p(2,:);
%          m = poly2mask(y,x, si(1), si(2));
%          m = imdilate(m, strel('square', 3)); % poly2 mask is very restrictive
%          n = size(obj.p,2);
%          %add lines
%          %for i = 1:n
%          %   m = impixelline(m, obj.p(:,i), obj.p(:,mod(i,n)+1), 1);
%          %end
%          x = round(x); y = round(y);
%          ids = and(and(x > 0, x <= si(1)), and(y > 0, y <= si(2)));
%          x = x(ids); y = y(ids);
%          ids = sub2ind(si, x,y); % add corners
%          m(ids) = 1;

         m = polygonToMask(obj.contour, 'size', si);
      end
      
      function n = npixel(obj, si)
         n = total(obj.mask(si));
      end
      
      function o = overlap(obj, roi)
         o = polygonArea(polygonIntersection(obj.contour), roi.toPolygon);
      end
      

      % exract roi from an array / image
      function [d, sh] = extractData(obj, d)
         %
         % [d, sh] = extractdata(obj, d)
         %
         % description:
         %     extracts data from bouding box
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
         pol = polygonShift(obj.contour, -shr);
         m = polygonToMask(pol, size(d));
         d = immask(d, m);
      end

      function shift(obj, sh)
         for i = 1:length(obj)
            obj(i).contour = polygonShift(obj(i).contour, sh);
         end
      end
      
      function obj = rescale(obj, scale)
         for i = 1:length(obj)
            obj(i).contour = polygonScale(obj(i).contour, scale);
         end
      end

      function obj= dilate(obj, bwidth)
         obj.p = polygonDilate(obj.contour, bwidth);
      end
      
      
      function plot(obj, varargin)
         n = length(obj);
         cc = colorcube(n);
         for i = 1:n
            polygonPlot(obj(i).contour, 'FaceColor', 'none', 'EdgeColor', cc(i,:), varargin{:})
         end
      end
      
  
   end
   
end