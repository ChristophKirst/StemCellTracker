classdef ROI < matlab.mixin.Copyable
   %
   % Shape class representing an abstract geometical shape in 2D or 3D eucleadian space
   %    
   
   
   methods
 
      function d = dim(obj)
         %
         % d = dim(obj)
         %
         % descritpion:
         %    dimension of the embedding space
        
         error('%s: dim not implemented!', class(obj));
      end
      
      function v = volume(obj, varargin)
         %
         % v = volume(obj)
         %
         % description
         %     volume of the shape in the embedding space
         
         error('%s: volume not implemented!', class(obj));
      end
      
      function b = boundingBox(obj, varargin)
         %
         % b = boundingBox(obj)
         %
         % description
         %     bounding box of the shape, returns a ROIRectangle
         
         error('%s: boundingBox not implemented!', class(obj));
      end

      function p = polygon(obj, varargin)
         %
         % p = polygon(obj, varargin)
         %
         % description
         %     polygonal representation / approximation of a 2D shape
         %     triangulation for a 3D shape
         %     typically the number of reference points needs to be specified
         %
         
         error('%s: polygon not implemented!', class(obj));
      end
      
      function m = mask(obj, varargin)
         %
         % m = mask(obj, si)
         %
         % description
         %    returns mask of ones of the regions pixels with the image of size si
         %

         error('%s: mask not implemented!', class(obj));
      end

      function n = nPixel(obj, varargin)
         %
         %  n = nPixel(obj, varargin)
         %
         % description
         %    returns number of pixel in the mask
         %

         n = total(obj.mask); % simple but slow / mem inensive
      end
      
      function plotRescaled(obj, scale, varargin)
         roi = obj.copy();
         roi.rescale(scale);
         plot(roi)
      end

   end
   
end