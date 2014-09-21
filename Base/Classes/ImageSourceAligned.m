classdef ImageSourceAligned < Alignment & ImageSource
   %
   % ImageSourceAligned class represents Image data coposed of aligned images
   % 

   methods   
      function obj = ImageSourceAligned(varargin)  % basic constructor
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
         obj.iinfo = obj.iinfo.fromData(data);
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

      function obj = setColor(obj, col)
         obj.initializeInfo();
         if iscell(col)
            obj.iinfo.icolor = col;
         else
            obj.iinfo.icolor = num2cell(col);
         end
      end
   
      function obj = setName(obj, nm)
         obj.initializeInfo();
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
      

      
      
      
      
      
      
      
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% alignment routines

      function sh = imageShifts(obj)  % image shifts from pairwise shifts
         %
         % shifts = imageShifts()
         %
         % description
         %      image shifts from pairwise shifts

         sh = obj.ialignment.imageShifts;
         
         if ~isempty(obj.itileformat)
            per = imuvwformat2permute('uvw', obj.itileformat);
            ts = obj.itileshape; ts = ts(per);
            sh = reshape(sh, ts);
            %sh = imuvwpermute(sh, obj.itileformat, 'uvw');
         end
      end

      function obj = alignPairsFromShifts(obj, ishifts)
         %
         % obj = alignPairsFromShifts(obj, ishifts)
         %
         % description:
         %    sets the pairwise shifts form image shifts
         %
         
         obj.ialignment.alignPairsFromShifts(ishifts);   
      end
      
      function obj = absoluteShiftsAndSize(obj)
         %
         % obj = absoluteShiftsAndSize(obj)
         %
         % description:
         %    calculates absolute size and shifts
         %
         % See also: absoluteShiftsAndSize
         
         obj.ialignment.absoluteShiftsAndSize(obj);
      end


      function obj = optimizePairwiseShifts(obj)
         %
         % obj = optimizePairwiseShifts(obj)
         %
         % description:
         %    globally optimizes pairwise shifts
         
         obj.ialignment.optimizePairwiseShifts;
      end
      
      function obj = makeShiftsConsistent(obj)
         %
         % obj = makeShiftsConsistent(obj)
         %
         % description:
         %    makes shifts mutually consistent (i.e. paths in the grid commute)
         
         obj.ialignment.makeShiftsConsistent;
      end
      
      
            
      function q = overlapQuality(obj, varargin)
         %
         % obj = overlapQuality(obj)
         %
         % description:
         %    calculates operlap quality of the images
         %
         % See also: overlapQuality, overlapStatisticsImagePair
         
         obj.ialignment.overlapQuality(obj, varargin{:});
         q = [obj.ialignment.ipairs.iquality];
      end

      function comp = connectComponents(obj, varargin)
         comp = obj.ialignments.connectedAlignments(varargin{:});
      end
      
      
      
      function obj = alignPairs(obj, varargin)
         %
         % alignPairs(obj, varargin)
         %
         % descritpion:
         %   alignes the individual paris of images
         %
         % input:
         %   param  parameter as for alignImagePair
         %
         % See also: alignImagePair
         
         obj.ialignment.alignPairs(obj, varargin{:});
      end
  
      
      function obj = align(obj, varargin)
         %
         % obj = align(obj, varargin)
         %
         % description:
         %    aligns images and sets new shifts
         %
         % See also: alignImages
         
         obj.ialignment.align(obj, varargin{:});
      end

      
      function st = stitch(obj, varargin)
         %
         % st = stitch(obj, source, param)
         %
         % description
         %     stitches images using the alignment information and source
         
         st = obj.ialignment.stitch(obj, varargin{:});
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