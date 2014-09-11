classdef ImageKeys < matlab.mixin.Copyable
   % keys to image subsets in an ImageSource
   
   properties
      ikeys = containers.Map;  % the map -> label and tag specs
   end

   methods
      function obj = ImageKeys(varargin)  % basic constructor
         %
         % ImageKeys()
         % ImageKeys(keys)
         %
         if nargin == 1 
            if isa(varargin{1}, 'ImageKeys') %% copy constructor
               obj = copy(varargin{1});
            elseif isa(varargin{1}, 'containers.Map')
               obj.icoordnames = varargin{1};
            else
               error('%s: not valid arguments for constructor', class(obj));
            end
         else
            error('%s: not valid arguments for constructor', class(obj));
         end  
      end
      
      function obj = setCoordinateNames(obj, names, coords)
         obj.icoordnames = containers.Map(names, coords);
      end
      
      function obj = addCoordinateName(obj, name, key)
         obj.icoordnames(name) = key;
      end

      function obj = removeCoordinateName(obj, name)
         obj.ikeys.remove(name);
      end
      
      function ids = indices(obj, name, iformat)
         ids = obj.icoordnames(name).indices(iformat);
      end

      function img = data(obj, varargin)
         img = obj.icoordnames(name).img(varargin{:});
      end
   end
 
   
end