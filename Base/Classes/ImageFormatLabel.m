classdef ImageFormatLabel < matlab.mixin.Copyable
  
   properties
      ilabel = containers.Map;  % the map label -> array of ImageFormatIndex classes
   end

   methods
      function obj = ImageFormatLabel(varargin)  % basic constructor
         %
         % ImageFormatLabel()
         % ImageFormatLabel(map)
         %
         if nargin == 0
            return
         elseif nargin == 1 
            if isa(varargin{1}, 'ImageFormatLabel') %% copy constructor
               obj = copy(varargin{1});
            elseif isa(varargin{1}, 'containers.Map')
               obj.ilabel = varargin{1};
            else
               error('%s: not valid arguments for constructor', class(obj));
            end
         else
            for l = 1:2:nargin
               if ~ischar(varargin{l})
                  error('%s: expects char array for label.', class(obj));
               end
               
               if iscell(varargin{l+1})
                  v2 = ImageFormatIndex(varargin{l+1});
               else
                  v2 = varargin{l+1};
               end
                  
               obj.addLabel(varargin{l}, v2);
            end
         end  
      end
      
      function obj = setLabels(obj, labels, specs)
         obj.ilabel = containers.Map(labels, specs);
      end
      
      function obj = addLabel(obj, label, spec)
         obj.ilabel(label) = spec;
      end

      function obj = removeLabel(obj, name)
         obj.ilabel.remove(name);
      end
      
      function ids = indices(obj, label, iformat)
         ids = obj.ilabel(label).indices(iformat);
      end

      function img = data(obj, label, varargin)
         img = obj.ilabel(label).data(varargin{:});
      end
   end
 
   
end