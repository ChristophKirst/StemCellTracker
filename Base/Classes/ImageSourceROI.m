classdef ImageSourceROI < ImageSource
   %
   % ImageSourceROI class represents a regoin of interest in a given image source
   % 

   properties 
      isource = [];    % image source
      iroi    = ROI;   % region of interest class
   end

   methods
      function obj = ImageSourceROI(varargin) % constructor
         %
         % ImageSourceTiled()
         % ImageSourceTiled(...,fieldname, fieldvalue,...)
         %
         
         if nargin == 0
            return
         elseif nargin == 1
            if isa(varargin{1}, 'ImageSourceROI') %% copy constructor
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

   end
   
   
end