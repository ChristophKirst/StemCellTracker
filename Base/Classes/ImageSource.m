classdef ImageSource < matlab.mixin.Copyable
   %
   % ImageSource class represents abstract Image data 
   % 
   % decription:
   %     this class can be used to represent an Image independent of its actual source
   %     use obj.data to return the actual image data
   %
   % required functions: size, format, class, color, tagmap
   %
   
   % Note: for translation to python, the cache structure should be via a global imge cache class
   %       that links ImageSource classes to the cached image
   %
   properties   
      iinfo  = [];          % info about image with properties/fields: .isize, .iformat, .icolors, .iclass, .icoords
      
      icache = false;       % (optional) weather to cache the data or not
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
   
            % obj.initialize();           
         end
      end
      
      function initialize(obj, varargin)
         obj.iinfo = obj.getInfo();
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % base methods
      %
      % usually not overwriten by super class
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
         
      function s = size(obj, varargin)
         if isempty(obj.iinfo.isize)
            s = obj.getSize(varargin{:});
            obj.iinfo.isize = s;
         else 
            s = obj.iinfo.isize;
         end
      end
      
      function d = ndims(obj)
         d = length(obj.size);
      end
      
      function d = sdims(obj)
         d = imsdims(obj.format);
      end
  
      function f = format(obj, varargin)
         if isempty(obj.iinfo.iformat)
            f = obj.getFormat(varargin{:});
            obj.iinfo.iformat = f;
         else 
            f = obj.iinfo.iformat;
         end
      end
   
      function c = class(obj)
         if isempty(obj.iinfo)
            obj.iinfo = obj.getInfo();
         end
         c = obj.iinfo.iclass;
      end
       
      function c = color(obj)
         if isempty(obj.iinfo)
            obj.iinfo = obj.getInfo();
         end
         c = obj.iinfo.icolor;
      end
      
      function s = scale(obj)         
         if isempty(obj.iinfo)
            obj.iinfo = obj.getInfo();
         end
         s = obj.iinfo.iscale;
      end

      function l = label(obj, varargin)
         if isempty(obj.iinfo)
            obj.iinfo = obj.getInfo();
         end
         l = obj.iinfo.ilabel;
      end 

      function i = info(obj)
         if isempty(obj.iinfo)
            obj.iinfo = obj.getInfo();
         end
         i = obj.iinfo;
      end
      
      
      % following base methods do not test for caching out of lazyness!
      function p = sizeP(obj)
         p = obj.iinfo.sizeP;
      end
      function q = sizeQ(obj)
         q = obj.iinfo.sizeQ;
      end
      function l = sizeL(obj)
         l = obj.iinfo.sizeL;
      end
      function c = sizeC(obj)
         c = obj.iinfo.sizeC;
      end
      function t = sizeT(obj)
         t = obj.iinfo.sizeT;
      end
      
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % access methods methods 
      % to obtain non-cached cached/preset info
      %
      % usually overwritten by sepcific super class
      
      function d = getData(obj, varargin)  % obtain the image data
         d = obj.idata;     % trivial here ->  to be implemented depending on ImageSource subclass
      end
      
      function obj = setData(obj, varargin)  % set the image data
         if nargin > 1
           obj.idata = varargin{1};
         end
      end
 
      function s = getSize(obj, varargin)
         s = size(obj.data);
      end
     
      function obj = setSize(obj, varargin)
         if nargin > 1
           obj.iinfo.isize = varargin{1};
         end
      end
     
     
      function f = getFormat(obj, varargin)
         f = imsize2format(obj.size);
      end
     
      function obj = setFormat(obj, varargin)
         if nargin > 1
           obj.iinfo.iformat = varargin{1};
         end
      end
     
      function c = getColor(varargin)
         c = {'gray'};
      end
     
      function obj = setColor(obj, varargin)
         if nargin > 1
            if iscell(varargin{1})
               obj.iinfo.icolor = varargin{1};
            else
               obj.iinfo.icolor = varargin(1);
            end
         end
      end
     
     
      function c = getClass(obj, varargin)
         c = class(obj.data);
      end
     
      function obj = setClass(obj, varargin)
         if nargin > 1
           obj.iinfo.iclass = varargin{1};
         end
      end
      
      function n = getName(obj, varargin) %#ok<INUSD>
         n = '';
      end
     
      function obj = setName(obj, varargin)
         obj.iinfo.iname = varargin{1};
      end

      function obj = setLabel(obj, varargin)
         obj.iinfo.ilabel = varargin{1};
      end
     
     
      function i = getInfo(obj, varargin)
         % get all info at once
         d = obj.data;
         
         i = ImageInfo();
         i.isize   = size(d);
         i.iformat = imsize2format(i.isize);
         i.iclass  = class(d);
         i.icolor  = obj.getColor();
         i.pqlctsizeFromFormatAndSize;
      end
     
      function obj = setInfo(obj, varargin)
         if nargin > 1
           obj.iinfo = varargin{1};
         end
      end

      % label extraction
      function d = extract(obj, label)
         if ~isempty(obj.iinfo.ilabel)
            d = obj.iinfo.ilabel.data(label, obj.data, obj.format);
         else
            warning('%s: cannot find label %s, return full image data!', class(obj), label);
            d = obj.data;
         end
      end

      % information / visulaization 
      function info = infoString(obj)
         cls = class(obj);
         info = ['ImageSource: ',  cls(12:end)];
         info = [info, '\nsize:   ', var2char(obj.size)];
         info = [info, '\nformat: ', var2char(obj.format)];
         info = [info, '\nclass:  ', var2char(obj.class)];
         info = [info, '\ncolor:  ', var2char(obj.color)];

%          if ~isempty(obj.label)
%             info = [info '\ntagmap:  '];
%             keys = obj.ilabel.keys;
%             vals = obj.ilabel.values;
%             for i = 1:obj.tagmap.length
%                info = [info, var2char(keys{i}), ' -> ' var2char(vals{i})]; %#ok<AGROW>
%             end
%          end
      end

      function print(obj)
         fprintf([obj.infoString, '\n']);
      end
      
      
      function plot(obj)
         if obj.sizeC == 1
            imcolormap(obj.color{1});
            implot(obj.data);
         elseif obj.sizeC == 3
            if isequal(obj.color, {'r', 'g', 'b'})
               implot(obj.data)
            else
               % TODO: handle special color settings here
               implot(obj.data);
            end
         else
            % TODO: handle non-standard channels here
            implot(obj.data);
         end
         
         if isentry(obj.iinfo, 'iname')
            set(gcf, 'Name', var2char(obj.iinfo.iname));
         end 
      end

   end
   
end