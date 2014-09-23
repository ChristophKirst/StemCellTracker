classdef ImageSource < matlab.mixin.Copyable
   %
   % ImageSource class represents abstract Image data 
   % 
   % decription:
   %     this class can be used to represent an Image independent of its actual source
   %     use obj.data to return the actual image data
   %
   % required functions: datasize, dataformat, color, data
   %
   
   % Note: for translation to python, the cache structure should be via a global imge cache class
   %       that links ImageSource classes to the cached image
   %
   properties   
      iinfo  = ImageInfo;   % info about image with properties/fields: .isize, .iformat, .icolors, .iclass, .icoords
      
      icache = true;        % (optional) weather to cache the data or not
      idata  = [];          % (optional) (cached) image data
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
               obj = obj.fromData(varargin{1});
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

      function obj = fromData(obj, data)
         obj.iinfo = ImageInfo();
         obj.iinfo.fromData(data);
         obj.idata = data; 
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % base methods
      %      

      function s = datasize(obj)
         s = obj.iinfo.idatasize;
      end
 
      function f = dataformat(obj)
         f = obj.iinfo.idataformat;
      end
       
      function s = rawsize(obj)
         s = obj.iinfo.irawsize;
      end
      
      function f = rawformat(obj)
         f = obj.iinfo.irawformat;
      end
      
         
      function d = subdata(obj, varargin)
         %
         % d = subdata(obj, datasepc)
         %
         % description:
         %    extract a subset of the data give the data specifications datasepc
         %
         
         param = parseParameter(varargin);
         dfrmt = num2cell(obj.dataformat);
         
         %find ids in param that match dfrmt
         fns = fieldnames(param);
         [ids, pos] = ismember(fns, dfrmt);
         pos = pos(ids);
         ids = find(ids);
            
         d = obj.data;
         asgn = repmat({':'}, ndims(d));
         for i = 1:length(ids)
            asgn{pos(i)} = param.(fns{ids(i)});
         end
         
         d = d(asgn{:});
      end
      
      function d = extractdata(obj, roi)
         %
         % d = subdata(obj, datasepc)
         %
         % description:
         %    extract a subset of the data given the spatial roi 
         %
         
         d = roi.extractdata(obj.data); 
      end
      
      function d = data(obj, varargin)
         if obj.icache % basic caching
            if ~isempty(obj.idata)
               d = obj.idata;
            else
               %get data the data
               d = obj.getData(varargin{:});
               
               % cache it
               obj.idata = d;
            end
         else
            d = obj.getData(varargin{:});
         end
      end
      
      function d = rawdata(obj, varargin)
         d = obj.getRawData();
      end
      
      
      function c = celldata(obj, varargin)
         c = {obj.data};  % trivial here
      end
      
      function cs = cellsize(obj, varargin)
         cs = obj.iinfo.icellsize;
      end

      function cf = cellformat(obj, varargin)
         cf = obj.iinfo.icellformat;
      end
      
      function n = ncells(obj, varargin)
         n = prod(obj.cellsize(varargin{:}));
      end
      
      
      function d = ndatadims(obj)
         % full number of image dimensinos
         d = length(obj.datasize);
      end
      
      function d = nsdatadims(obj)
         % spatial image dimensions
         d = imsdims(obj.dataformat);
      end

      function c = dataclass(obj)  % we dont use class here as this would overwrte the class routine
         c = obj.iinfo.idataclass;
      end
       
      function c = color(obj)
         c = obj.iinfo.icolor;
      end
      
      function s = scale(obj)         
         s = obj.iinfo.iscale;
      end
      
      function l = name(obj)
         l = obj.iinfo.iname;
      end 

      function i = info(obj)
         i = obj.iinfo;
      end
      
      function s =datasizePQLCT(obj)
         s = obj.iinfo.idatasizePQLCT;
      end
      

      function p = datasizeP(obj)
         p = obj.iinfo.datasizeP;
      end
      function q = datasizeQ(obj)
         q = obj.iinfo.datasizeQ;
      end
      function l = datasizeL(obj)
         l = obj.iinfo.datasizeL;
      end
      function c = datasizeC(obj)
         c = obj.iinfo.datasizeC;
      end
      function t = datasizeT(obj)
         t = obj.iinfo.datasizeT;
      end

      
      function i = dataformatpos(obj, frmt)       
         i = obj.iinfo.dataformatpos(frmt);
      end
      
      function i = rawformatpos(obj, frmt)
         i = obj.iinforawformatpos(frmt);
      end
      
      function i = cellformatpos(obj, frmt)
         i = obj.iinfo.cellformatpos(frmt);
      end


      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % access methods 
      % to obtain non-cached/preset info
      %
      % usually overwritten by sepcific super class

      function d = getData(obj, varargin)  % obtain the image data
         %
         % d = getData(obj, varargin)
         %
         % description:
         %     obtain image data in the format given by obj.dataformat
         
         d = obj.idata;     % trivial here ->  to be implemented depending on ImageSource superclass
      end
      
      function d = getRawData(obj, varargin)
         %
         % d = getRawData(obj, varargin)
         %
         % description:
         %     obtain raw image data in the format given by obj.rawformat
         
         d = obj.iinfo.data2raw(obj.idata);     % trivial here ->  to be implemented depending on ImageSource superclass
      end
      
      function obj = setData(obj, d)  % set the image data
         obj.idata = d;
      end

      function obj = setDataFormat(obj, dfrmt)         
         obj.clearCache();
         obj.iinfo.idataformat = dfrmt;
      end
       
      function obj = setColor(obj, col)
         if iscell(col)
            obj.iinfo.icolor = col;
         else
            obj.iinfo.icolor = num2cell(col);
         end
      end
   
      function obj = setName(obj, nm)
         obj.iinfo.iname = nm;
      end
      
      
      function obj = setInfo(obj, info)
         obj.iinfo = info;
      end

 
      function obj = clearCache(obj)
         obj.idata = [];
      end
      
      function obj = setCache(obj, c)
         obj.icache = c;
         obj.clearCache();
      end
      


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
         istr = [istr, '\nsize:   ', var2char(obj.datasize)];
         istr = [istr, '\nformat: ', var2char(obj.dataformat)];
         istr = [istr, '\nclass:  ', var2char(obj.dataclass)];
         istr = [istr, '\ncolor:  ', var2char(obj.color)];
      end

      function print(obj)
         fprintf([obj.infoString, '\n']);
      end

      function plot(obj)
         implotis(obj);
      end
      
   end
   
end