classdef ImageScale < matlab.mixin.Copyable 
   %
   % ImageScale class
   %
   % description:
   %     class representing spatial information of the image
   %
   % See also: ImageInfo

   properties
      scale   = [1; 1];  % spatial extend of one pixel in units given by unit
      units   = 'um';   % unit
   end

   methods
      function obj = ImageScale(varargin)  % basic constructor
         %
         % ImageScale()
         % ImageScale(scale)
         % ImageScale(...,fieldname, fieldvalue,...)
         %
         if nargin == 1 
            if isa(varargin{1}, 'ImageScale') %% copy constructor
               obj = copy(varargin{1});
            elseif isnumeric(varargin{1})
               obj.scale = varargin{1};
            else
               error('%s: not valid arguments for constructor', class(obj));
            end
         else
            for i = 1:2:nargin % constructor from arguments
               if ~ischar(varargin{i})
                  error('%s: invalid constructor input, expects char at position %g',class(obj), i);
               end
               if isprop(obj, varargin{i})
                  obj.(varargin{i}) = varargin{i+1};
               else
                  warning('%s: unknown property name: %s ', class(obj), varargin{i})
               end
            end
   
            obj.initialize();           
         end
      end
      
      function intialize(obj)
         obj.scale = obj.scale(:);
      end
      
      
      function s = pixel2space(obj, p)
         s = bsxfun(@times, p, obj.scale);
      end

      function p = space2pixel(s)
         p = bsxfun(@rdivide, s, obj.scale);
      end
      
      
   end
