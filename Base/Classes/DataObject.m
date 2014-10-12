classdef DataObject < Object
%
% DataObject class extends the Object class with a data and segmentation property 
% for object specific data such as fluorescence measurement, or segmentation results
%
% See also: Object, Segment, Frame, Trajectory, TimeSeries
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
               obj.fromObject(varargin{1});
            end
         else
            obj.fromParameter(varargin);
         end
      end
      
      function obj = fromParameter(obj, varargin)
         obj = classFromParameter(obj, [], varargin);
      end
      
      function obj = fromObject(obj, o) 
         oldobj = o;
         nobjs = length(oldobj);
         obj(nobjs) = DataObject();
         props = properties(oldobj);
         for k = 1:nobjs
            for i = 1:length(props)
               obj(k).(props{i}) = oldobj(k).(props{i});
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
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % basics

      % get statistics names
      function fnames = dataFields(obj)
         %
         % fnames = fields(obj);
         %
         % description:
         %     returns field names of the data entries

         if isempty(obj(1).data)
            fnames = {};
         else
            fnames = fieldnames(obj(1).data);
         end
      end
      
      
      function val = dataValues(obj, sname)
         %
         % val = obj.dataValue(sname)
         %
         % description:
         %     returns the value of the field sname
            
            val = [obj.data];
            val = [val.(sname)]; 
      end
      
      
      function cname = channelName(~, ch, fname)
      %
      % val = channelname(ch, fname)
      %
      % description:
      %     returns the value of the field ['ch' num2str(ch) '_' fname]
      
         if nargin < 3
            fname = '';
         end
         if isempty(fname)
            cname = ['ch_', num2str(ch)];
         else
            cname = ['ch_', num2str(ch) '_' fname];
         end
      end
      
      
      % easier data access
      function val = channelData(obj, varargin)
      %
      % val = channel(ch, fname)
      %
      % description:
      %     returns the value of the field ['ch' num2str(ch) '_' fname]
         
         cname = obj.channelName(varargin{:});
         
         val = [obj.data];
         val = [val.(cname)];
      end
      
      % set channel data
      function obj = setChannelData(obj, ch, varargin)
      %
      % addChannelData(obj, ch, stats)
      % addChannelData(obj, ch, values)
      % addChannelData(obj, ch, fname, values)
      % addChannelData(obj, ch, stats, fnames)
      %
      % description:
      %     adds/overwrites the entries in stats with prefix ch to the obj.data structure

         if isstruct(varargin{1})
            stats = varargin{1};
            
            if nargin < 4
               fnames = fieldnames(stats);
            else
               fnames = varargin{2};
            end
            if ischar(fnames)
               fnames = {fnames};
            end
            
            if isempty(fnames)
               return
            end

            for i = length(fnames):-1:1
               newfnames{i} = ['ch_', num2str(ch) '_' fnames{i}];
            end

            stats = rmfield(stats, setdiff(fieldnames(stats), fnames));
            stats = renamestruct(stats, fnames, newfnames);
            obj = obj.setData(stats);
            
         elseif ischar(varargin{1}) 
            
            fname = varargin{1};
            vals  = varargin{2};
            cname = obj.channelName(ch, fname);
            obj.setData(cname, vals);

         else
            
            vals  = varargin{1};
            cname = obj.channelName(ch);
            obj.setData(cname, vals);

         end
 
      end
      
      
      %add statistics
      function obj = setData(obj, varargin)
         %
         % setData(obj, stats)
         % setData(obj, fname, data)
         %
         % description:
         %     adds the entries in stats to the obj.data structure, or adds a field fname with data to obj.data
         
         if nargin <= 2
            if nargin == 0 || ~isstruct(varargin{1})
               error('%s: addData expects struct as input!', class(obj))
            end
            
            stats = varargin{1};

            fnames = obj(1).dataFields();   % assume identical field names !
            snames = fieldnames(stats);
              
            if ~isempty(fnames)
               dstats = [obj.data];
               dstats = rmfield(dstats, intersect(fnames, snames));
            else
               dstats(length(obj)) = struct();
            end

            for i = 1:length(snames)
               dd = num2cell([stats.(snames{i})]);
               [dstats.(snames{i})] = dd{:};
            end
            dstats = num2cell(dstats);
            [obj.data] = dstats{:};
         else
            if ~ischar(varargin{1})
               error('%s: addData expects data name as first argument!', class(obj))
            end
            
            fname = varargin{1};
            
            dstats = [obj.data];
            d = num2cell(varargin{2});
            [dstats.(fname)] = d{:};
            
            dstats = num2cell(dstats);
            [obj.data] = dstats{:};
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
      
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % specialized data access functions
      
      function val = dapi(obj)
         %
         % val = dapi(obj)
         %
         % description:
         %     returns the value of the field 'dapi'
            val = [obj.data];
            val = [val.('ch_dapi')]; 
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

      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % visualization
      
      function plot(obj)
         implot(obj.labeledImage)
      end
      

   end % methods
end % classdef
      
      