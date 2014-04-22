classdef DataObject < Object
%
% DataObject class extends the Object class with a data and segmentation property 
% for object specific data such as fluorescence measurement, or segmentation results
%
% See also: Cell, Frame, Trajectory
   properties
      data = [];
      segment = [];
   end
   
   methods
      function obj = DataObject(varargin)  % simple constructor
      %
      % DataObject()
      % DataObject(...,fieldname, fieldvalue,...)
      %
         if nargin == 1
            if isa(varargin{1}, 'DataObject') %% copy constructor
               obj = copy(varargin{1});
            elseif isa(varargin{1}, 'Object')
               oldobj = varargin{1};
               nobjs = length(oldobj);
               obj(nobjs) = DataObject();
               props = properties(oldobj);
               for k = 1:nobjs
                  for i = 1:length(props)
                     obj(k).(props{i}) = oldobj(k).(props{i});
                  end
               end
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
      end

      function newobj = copy(obj)
         %
         % c = copy(obj)
         %
         % description:
         %     deep copy of the object
         
         nobjs = length(obj);
         newobj(nobjs) = DataObject();
         props = properties(newobj);
         for k = 1:nobjs
            for i = 1:length(props)
               newobj(k).(props{i}) = obj(k).(props{i});
            end
         end
      end
      
      
      % get statistics names
      function fnames = dataFields(obj)
         %
         % fnames = fields(obj);
         %
         % description:
         %     returns field names of the data entries

         fnames = fieldnames(obj(1).data);    
      end
      
      
      function val = stat(obj, sname)
         %
         % val = stat(sname)
         %
         % description:
         %     returns the value of the field ['ch' num2str(ch) '_' fname]
            
            val = [obj.data];
            val = [val.(sname)]; 
      end
      
      
      % easier data access
      function val = channel(obj, ch, fname)
         %
         % val = channel(ch, fname)
         %
         % description:
         %     returns the value of the field ['ch' num2str(ch) '_' fname]
            
            fname = ['ch' num2str(ch) '_' fname];
            val = [obj.data];
            val = [val.(fname)]; 
      end
      
      %some specialized data access function 
      function val = dapi(obj)
         %
         % val = dapi(obj)
         %
         % description:
         %     returns the value of the field 'dapi'
            val = [obj.data];
            val = [val.('dapi')]; 
      end
      
      %add more here if needed...


      % objects segmentation
      function seg = nucleus(obj)
         seg = arrayfun(@first, [obj.segment], 'UniformOutput', false);
         seg = [seg{:}];
      end
      
      function seg = cytoplasm(obj)
         if size(obj(1).segment,1) > 1
            seg = arrayfun(@last, [obj.segment]);
         else
            seg = [];
         end
      end

      function imglab = labeledImage(obj)
         shape = [obj.segment];
         [idx, isize] = shape.pixelIdxList();
         imglab = zeros(isize);
         for i = 1:length(idx)
            imglab(idx{i}) = obj(i).id;
         end
      end
      

   end % methods
end % classdef
      
      