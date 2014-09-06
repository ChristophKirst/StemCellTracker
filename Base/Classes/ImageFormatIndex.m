classdef ImageFormatIndex < matlab.mixin.Copyable
  
   properties
      icoordinate = '';   % the image format coordinate name (e.g. 'c')
      iindex      = [];   % the index range, including ':' for all
   end
   
   methods
      function obj = ImageFormatIndex(varargin)  % basic constructor
         %
         % ImageFormatIndex(iformat, iindex)
         %
         if nargin == 0
            return
         elseif nargin == 1 
            if isa(varargin{1}, 'ImageFormatIndex') %% copy constructor
               obj = copy(varargin{1});
            elseif isa(varargin{1}, 'cell')
               n = ceil(length(varargin{1})/2);
               obj(n) = ImageFormatIndex;
               obj.fromCell(varargin{1});
            else
               error('%s: not valid arguments for constructor', class(obj));
            end
         else
            n = ceil(length(varargin)/2);
            obj(n) = ImageFormatIndex;
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
            error('%s: inconsistent arguments for fromKeys', class(obj));
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
   end
end