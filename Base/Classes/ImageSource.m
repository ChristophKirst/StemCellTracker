classdef ImageSource < matlab.mixin.Copyable
   %
   % ImageSource class represents abstract Image data 
   % 
   % decription:
   %     this class can be used to represent an Image independent of its actual source
   %     use obj.data to return the actual image data
   %
   % required functions: size, format, class, color, label
   %
   
   % Note: for translation to python, the cache structure should be via a global imge cache class
   %       that links ImageSource classes to the cached image
   %
   properties   
      iinfo  = [];          % info about image with properties/fields: .isize, .iformat, .icolors, .iclass, .icoords
      
      icache = true;        % (optional) weather to cache the data or not
      idata  = [];          % (optional) cached image data
   end
        
   methods   
      function obj = ImageSource(varargin)  % basic constructor
         %
         % ImageSource()
         % ImageSource(...,fieldname, fieldvalue,...)
         %
         if nargin == 1 
            if isa(varargin{1}, 'ImageSource') %% copy constructor
               obj = copy(varargin{1});
            elseif isnumeric(varargin{1})
               obj.idata = varargin{1};
               obj.initialize;
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
   
            %obj.initialize();           
         end
      end
      
      function initialize(obj, varargin)
         obj.iinfo = obj.getInfo();
      end
      
      function initializeInfo(obj)
         if isempty(obj.iinfo)
            obj.iinfo = obj.getInfo();
         end
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % base methods
      %

      function d = data(obj, varargin)
         if obj.icache % basic caching
            if ~isempty(obj.idata)
               d = obj.idata;
            else
               d = obj.getData(varargin{:});
               obj.idata = d;
            end
         else
            d = obj.getData(varargin{:});
         end
      end
      
      function c = celldata(obj, varargin)
         c = {obj.data};  % trivial here
      end
      
      function cs = cellsize(obj, varargin)
         cs = obj.iinfo.icellsize;
      end
         
      function s = size(obj)
         %obj.initializeInfo;
         s = obj.iinfo.isize;
      end
      
      function d = ndims(obj)
         d = length(obj.size);
      end
      
      function d = sdims(obj)
         d = imsdims(obj.format);
      end
  
      function f = format(obj)
         %obj.initializeInfo;
         f = obj.iinfo.iformat;
      end
   
      function c = class(obj)
         %obj.initializeInfo;
         c = obj.iinfo.iclass;
      end
       
      function c = color(obj)
         %obj.initializeInfo;
         c = obj.iinfo.icolor;
      end
      
      function s = scale(obj)         
         %obj.initializeInfo;
         s = obj.iinfo.iscale;
      end

      function l = keys(obj)
         %obj.initializeInfo;
         l = obj.iinfo.ikeys;
      end 
      
      function l = name(obj)
         %obj.initializeInfo;
         l = obj.iinfo.iname;
      end 

      function i = info(obj)
         %obj.initializeInfo;
         i = obj.iinfo;
      end
      
      

      function p = sizeP(obj)
         %obj.initializeInfo;
         p = obj.iinfo.sizeP;
      end
      function q = sizeQ(obj)
         %obj.initializeInfo;
         q = obj.iinfo.sizeQ;
      end
      function l = sizeL(obj)
         %obj.initializeInfo;
         l = obj.iinfo.sizeL;
      end
      function c = sizeC(obj)
         %obj.initializeInfo;
         c = obj.iinfo.sizeC;
      end
      function t = sizeT(obj)
         %obj.initializeInfo;
         t = obj.iinfo.sizeT;
      end

      
      % cell sizes for certain super classes
%       function u = sizeU(~)
%          u = 0;
%       end   
%       function v = sizeV(~)
%          v = 0;
%       end   
%       function w = sizeW(~)
%          w = 0;
%       end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % access methods 
      % to obtain non-cached/preset info
      %
      % usually overwritten by sepcific super class
      
      function d = getData(obj, varargin)  % obtain the image data
         d = obj.idata;     % trivial here ->  to be implemented depending on ImageSource superclass
      end
      function obj = setData(obj, varargin)  % set the image data
         obj.idata = varargin{1};
      end
 
      
      function s = getSize(obj, varargin)
         obj.initializeInfo();
         s = obj.iinfo.isize;
      end
      function obj = setSize(obj, varargin)
         obj.initializeInfo();
         obj.iinfo.isize = varargin{1};
      end
      
      
      function f = getFormat(obj, varargin)
         obj.initializeInfo();
         f = obj.iinfo.iformat;
      end
      function obj = setFormat(obj, varargin)
         obj.initializeInfo();
         obj.iinfo.iformat = varargin{1};
      end
     
      
      function c = getColor(varargin)
         obj.initializeInfo();
         c = obj.iinfo.icolor;
      end
      function obj = setColor(obj, varargin)
         obj.initializeInfo();
         if iscell(varargin{1})
            obj.iinfo.icolor = varargin{1};
         else
            obj.iinfo.icolor = varargin(1);
         end
      end
     
      
      function c = getClass(obj, varargin)
         obj.initializeInfo();
         c = obj.iinfo.iclass;
      end
      function obj = setClass(obj, varargin)
         obj.initializeInfo();
         obj.iinfo.iclass = varargin{1};
      end
      
      
      function n = getName(obj, varargin)
         obj.initializeInfo();
         n = obj.iinfo.iname;
      end
      function obj = setName(obj, varargin)
         obj.iinfo.iname = varargin{1};
      end
      
      
%       function k = getKeys(obj, varargin)
%          obj.initializeInfo();
%          k = obj.iinfo.ikeys;
%       end
%       function obj = setKeys(obj, varargin)
%          obj.initializeInfo();
%          obj.iinfo.ikeys = varargin{1};
%       end
%       
%       function k = getKey(obj, varargin)
%          obj.initializeInfo();
%          k = obj.iinfo.ikeys.getKey(varargin{:});
%       end
%       function obj = setKey(obj, varargin)
%          obj.initializeInfo();
%          obj.iinfo.ikeys.setKey(varargin{:});
%       end
     
      
      function i = getInfo(obj, varargin)
         % get all info at once         
         i = imdata2info(obj.data);
      end
      function obj = setInfo(obj, varargin)
         obj.iinfo = varargin{1};
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % data extraction using a key
%       
%       function d = extract(obj, label)
%          if ~isempty(obj.iinfo.ilabel)
%             d = obj.iinfo.ilabel.data(label, obj.data, obj.format);
%          else
%             warning('%s: cannot find label %s, returning full image data!', class(obj), label);
%             d = obj.data;
%          end
%       end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % information / visulaization 
      function istr = infoString(obj, varargin)
         if nargin > 1
            cls = varargin{1};
         else
            cls = '';
         end
         istr = ['ImageSource: ',  cls];
         if ~isempty(obj.name)
            istr = [istr, '\nname:   ', var2char(obj.name)];
         end
         istr = [istr, '\nsize:   ', var2char(obj.size)];
         istr = [istr, '\nformat: ', var2char(obj.format)];
         istr = [istr, '\nclass:  ', var2char(obj.class)];
         istr = [istr, '\ncolor:  ', var2char(obj.color)];
%         istr = [istr, '\n', obj.label.infoString];
      end

      function print(obj)
         fprintf([obj.infoString, '\n']);
      end
      
      
      function plot(obj)
         implotis(obj);
      end

   end
   
end