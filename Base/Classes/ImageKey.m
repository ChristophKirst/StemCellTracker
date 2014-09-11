classdef ImageKey < matlab.mixin.Copyable
  
   properties
      idatakey = [];   % data specs
      itagkey  = [];   % tag specs
   end
   
   methods
      function obj = ImageKey(varargin)  % basic constructor
         %
         % ImageLabelIndex(iformat, iindex)
         %
         if nargin == 0
            return
         elseif nargin == 1 
            if isa(varargin{1}, 'ImageLabelIndex') %% copy constructor
               obj = copy(varargin{1});
            elseif isa(varargin{1}, 'cell')
               n = ceil(length(varargin{1})/2);
               obj(n) = ImageLabelIndex;
               obj.fromCell(varargin{1});
            else
               error('%s: not valid arguments for constructor', class(obj));
            end
         else
            n = ceil(length(varargin)/2);
            obj(n) = ImageLabelIndex;
            obj.fromCell(varargin);
         end 
            
         obj.initialize();
         
      end
      
      function initialize(obj)
         for i = 1:length(obj)
            if ~isempty(obj(i).icoordinate)
               if ~ any(obj(i).icoordinate == 'xyzpqlct')
                  error('%s: format coordinate %s not in xyzpqlct', class(obj(i)), obj(i).icoordinate);
               end
            end
         end
      end
      
      function obj = fromCell(obj, ce)
         if mod(length(ce(:)),2) ~=0
            error('%s: not valid arguments for constructor', class(obj));
         end
  
         k = 1;
         for i = 1:2:length(ce)
            obj(k).icoordinate = ce{i};
            obj(k).iindex      = ce{i+1};
            k = k + 1;
         end  
      end
      
      
      function ids = indices(obj, iformat)
      % gives cell array such that img(ids{:}) extracts correct subset
         ids = repmat({':'},1, length(iformat));
         for i = 1:length(obj)
            pos = find(obj(i).icoordinate == iformat, 1);
            if ~isempty(pos)
               ids{pos} = obj(i).iindex;
            end
         end
      end

      function img = data(obj, varargin)
         if nargin < 2 || nargin > 3
            error('%s: inconsistent number of arguments', class(obj));
         elseif nargin == 2
            img = varargin{1};
            iformat = imformat(varargin{1});
         elseif nargin == 3
            img = varargin{1};
            iformat = varargin{2};
         end
         
         ids = obj.indices(iformat);
         img = img(ids{:});
      end
      
      
      function istr = infoString(obj)
         istr = '(';
         for i = 1:length(obj)
            istr = [istr, obj(i).icoordinate, '.', var2char(obj(i).iindex), ', ']; %#ok<AGROW>
         end
         istr = [istr(1:end-2), ')'];
      end
      
   end
end