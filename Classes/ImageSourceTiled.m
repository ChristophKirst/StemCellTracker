classdef ImageSourceTiled < ImageSource
   %
   % ImageSourceTiled class represents an image composed from a tiling
   % 

   properties 
      images = {};  % cell array of ImageSource classes specifiing the individual images
      shifts = {};  % relative shifts 
      alignparam  = [];  % parameter needed for stitching
      stichparam  = [];  % parameter needed for alignement
   end

   methods
      function obj = ImageSourceTiled(varargin) % constructor
         %
         % ImageSourceTiled()
         % ImageSourceTiled(...,fieldname, fieldvalue,...)
         %
         
         if nargin == 0
            return
         elseif nargin == 1
            if isa(varargin{1}, 'ImageSourceTiled') %% copy constructor
               obj = copy(varargin{1});
            elseif iscell(varargin{1})
               obj.images = varargin{1};
            else
               error('%s: invalid constructor input, expects cell array at position %g',class(obj), 1);
            end
         else
            for i = 1:2:nargin % constructor from arguments
               if ~ischar(varargin{i})
                  error('%s: invalid constructor input, expects char at position %g',class(obj), i);
               end
               if isprop(obj, lower(varargin{i}))
                  obj.(lower(varargin{i})) = varargin{i+1};
               else
                  warning('%s: unknown property name: %s ', class(obj), lower(varargin{i}))
               end
            end
         end

         %some automatic conversion form file names to ImageSourceFile objectgs etc
         
         
         if isempty(obj.images)
            error('%s: no images specified in constructor!', class(obj));
         end
         
         %if isempty(obj.shifts)
         %   obj.shifts = alignImages(obj.getTiledData(), obj.alignparam);
         %end

         %if isempty(obj.size)
         %   obj.size = obj.getSize();
         %end

         %if isempty(obj.format)
         %   obj.format = obj.getFormat();
         %end
      end

      function img = getRawData(obj, varargin)
         img = obj.getTiledData();
         img = stitchImages(img, obj.shifts, obj.stichparam);
      end
      
      function imgs = getTiledRawData(obj, varargin)
         imgs = cellfun(@(x) x.getRawData(varargin{:}), obj.images, 'UniformOutput', false); 
      end
  
      function imgs = getTiledData(obj, varargin)
         imgs = cellfun(@(x) x.getData(varargin{:}), obj.images, 'UniformOutput', false);
      end
     
      function isiz = getTiledImageSizes(obj, varargin)
         isiz = cellfun(@(x) x.size, obj.images, 'UniformOutput', false);
      end
         
      function siz = getSize(obj, varargin)   
         [~, siz] = absoluteShiftsAndSizes(obj.shifts, obj.getTiledImageSizes());
      end

      function ashifts = getAbsoluteShifts(obj, varargin)
         [ashifts, ~] = absoluteShiftsAndSizes(obj.shifts, obj.getTiledImageSizes());
      end
      
      function [ashifts, asiz] = getAbsoluteShiftsAndSizes(obj, varargin)
         [ashifts, asiz] = absoluteShiftsAndSizes(obj.shifts, obj.getTiledImageSizes());
      end
      
      function shifts = align(obj)
         shifts = alignImages(obj.getTiledData(), obj.alignparam);
         obj.shifts = shifts;
         obj.initialize();
      end

             
      %function it = getSubTiling(obj, varargin) % extrack a sub tiled image from this one
      %end

      function info = infoString(obj)
         info = infoString@ImageSource(obj);
         info = [info, '\ntiling: ', var2char(size(obj.images))];
         info = [info, '\nshifts: ', var2char(obj.shifts)];
      end

   end
   
   
end