classdef ImageLabels < matlab.mixin.Copyable
  
   properties
      ilabels = containers.Map;  % the map label -> array of ImageLabelIndex classes
   end

   methods
      function obj = ImageLabels(varargin)  % basic constructor
         %
         % ImageLabels()
         % ImageLabels(map)
         %
         if nargin == 0
            return
         elseif nargin == 1 
            if isa(varargin{1}, 'ImageLabels') %% copy constructor
               obj = copy(varargin{1});
            elseif isa(varargin{1}, 'containers.Map')
               obj.ilabels = varargin{1};
            else
               error('%s: not valid arguments for constructor', class(obj));
            end
         else
            for l = 1:2:nargin
               if ~ischar(varargin{l})
                  error('%s: expects char array for label.', class(obj));
               end
               
               if iscell(varargin{l+1})
                  v2 = ImageDataIndex(varargin{l+1});
               else
                  v2 = varargin{l+1};
               end
                  
               obj.addLabel(varargin{l}, v2);
            end
         end  
      end
      
      function obj = setLabels(obj, labels, specs)
         obj.ilabels = containers.Map(labels, specs);
      end
      
      function obj = addLabel(obj, label, spec)
         obj.ilabels(label) = spec;
      end

      function obj = removeLabel(obj, name)
         obj.ilabels.remove(name);
      end
      
      function ids = indices(obj, label, iformat)
         ids = obj.ilabels(label).indices(iformat);
      end

      function img = data(obj, label, varargin)
         img = obj.ilabels(label).data(varargin{:});
      end
      
      
      function istr = infoString(obj)
         istr = 'label:  ';
         if ~isempty(obj.ilabels)

            keys = obj.ilabels.keys;
            vals = obj.ilabels.values;
            for i = 1:length(keys)
               istr = [istr, var2char(keys{i}), ' -> ' vals{i}.infoString() ';  ']; %#ok<AGROW>
            end
         else
            istr = [istr 'none'];
         end
      end
      
      
   end
 
   
end