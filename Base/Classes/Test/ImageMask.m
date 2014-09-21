classdef ImageMask
   %
   % ImageMask class representing a mask of an image
   %
   
   properties
      itype  = 'mask' ;    % type of the mask: 'mask' = binary mask, 'shape' a geometric shape
      imask  = [];
   end
   
   methods
      function obj = ImageMask(varargin)  % basic constructor
         %
         % ImageMask()
         % ImageMask(...,fieldname, fieldvalue,...)
         %
         if nargin == 1 
            if isa(varargin{1}, 'ImageMask') %% copy constructor
               obj = copy(varargin{1});
            elseif isnumeric(varargin{1})
               obj.imask = varargin{1};
               obj.itype = 'maks';
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
         end
      end
      
      
      function m = mask(obj)
         switch obj.itype
            case 'mask'
               m = obj.imask;
            case 'shape'
               error('ImageMask: not iplemented yet!')
            otherwise
               error('ImageMask: unkonw masking type: %s!', var2char(mask))
         end
      end
   
   end 
end