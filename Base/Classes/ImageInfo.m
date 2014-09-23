classdef ImageInfo < matlab.mixin.Copyable
   
   % class that manages information about a single image in a file
   % for multiple series in a file choose ImageSourceList
   
   % notes:
   %    the info class is based on the info needed to read images via imread_bf / loci tools
   %    data read from files always linked to the pqlct dims to simplify array manipulations / reading etc
   %    series can be linked to cell or image dims
   %    tags can be linked to cell or image dims
   
   properties
      idatasize      = [];           % size of the image as returned by .data routine 
      idataformat    = '';           % format of the image as returned by .data routine
     
      idatasizePQLCT = [0,0,0,0,0];  % size of the image for full pqlct dimensions

      icellsize   = 1;               % size of the cell structure as returend when using .celldata routine
      icellformat = '';              % format of the cell structure (usually fixed subset of uvwrs)
     
      irawsize   = [];               % size of the raw data, [] = idatasize
      irawformat = '';               % format of the raw image data, when read from a file, '' = iformat
  
      idataclass = '';               % class of the image
            
      icolor = {};                   % colors for the channels
      ialpha = {};                   % (optional) alpha value or alpha mask for each channel (=1 if no specified) 
      
      iseries = 1;                   % (optional) series number or ids if image if from a larger file 
      iseriesformat = '';            % (optional) dimension the series is linked to 

      inimages = [];                 % (optional) number of total 2d images as counted via loci tools

      imetadata = [];                % (optional) meta data
      
      iscale = [];                   % (optional) spatial scale of image in pixel per spatial unit
      iunit = '';                    % (optional) spatial unit
      
      iname  = [];                   % (optional) name of the image 
   end
   
   methods
      
      
      function obj = ImageInfo(varargin) % constructor
         %
         % ImageInfo()
         % ImageInfo(data)
         % ImageInfo(...,fieldname, fieldvalue,...)
         %
         
         if nargin == 0
            return
         elseif nargin == 1
            if isa(varargin{1}, 'ImageInfo') %% copy constructor
               obj = copy(varargin{1});
            elseif isnumeric(varargin{1})
               obj.fromData(varargin{1});
            else
               error('%s: invalid constructor input, expects char at position %g',class(obj), 1);
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
         
         obj.icache = false; % for tagging simple caching is usaually not a good idea
      end
      
      
      function obj = fromData(obj, data)
         info = imdata2info(data);
         for n = properties(obj)'
            obj.(n{1}) = info.(n{1});
         end
      end

      function obj = fromFile(obj, filename)
         info = imread_bf_info(filename);
         for n = properties(obj)'
            obj.(n{1}) = info.(n{1});
         end
      end
      
      function obj = fromImfinfo(obj, info)
         info = iminfo2info(info);
         for n = properties(obj)'
            obj.(n{1}) = info.(n{1});
         end
      end
      
      function obj = pqlctsizeFromFormatAndSize(obj)
         frmt = 'pqlct';  
         pos = arrayfun(@(x) find(x == frmt,1), obj.idataformat);
         obj.idatasizePQLCT = ones(1, 5);
         obj.idatasizePQLCT(pos) = obj.idatasize;
      end
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% sizes
      
      function d = datadim(obj)
         d = length(obj.idataformat);
      end
      
      function d = celldim(obj)
         d = length(obj.icellformat);
      end
   
      
      function s = datasize(obj, varargin)
         if nargin == 1
            s = obj.idatasize;
         else
            s = obj.idatasize(obj.dataformatpos(varargin{1}));
         end
      end
      
             
      function si = dataformatsize(obj, frmt)
         %
         % si = dataformatsize(obj, frmt)
         %
         % description:
         %     returns the size of the dimensions specified by frmt
         %
         i  = obj.dataformatpos(impqlformat2format(frmt)); % for size inversion of the axes does not matter
         si = obj.idatasize;
         si = si(i);
      end
      
      function si = rawformatsize(obj, frmt)
         i = obj.rawformatpos(impqlformat2format(frmt));
         si = obj.irawsize;
         si = si(i);
      end

      
      function p = datasizeP(obj)
         p = obj.idatasizePQLCT(1);
      end
      function q = datasizeQ(obj)
         q = obj.idatasizePQLCT(2);
      end
      function l = datasizeL(obj)
         l = obj.idatasizePQLCT(3);
      end
      function c = datasizeC(obj)
         c = obj.idatasizePQLCT(4);
      end
      function t = datasizeT(obj)
         t = obj.idatasizePQLCT(5);
      end

     
      function obj = setDataSize(obj, newsize) %using this routine is usally not a good idea
         obj.idatasize = newsize;
         obj.pqlctsizeFromFormatAndSize();
      end
      

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% formatting / permuations

      function d = raw2data(obj, d, varargin)
         if nargin > 2
            rfrmt = varargin{1};
         else
            rfrmt = obj.irawformat;
         end
         
         if ~isempty(rfrmt)
            d = impqlpermute(d, rfrmt, obj.idataformat);
         end 
      end
      
      
      function d = data2raw(obj, d, varargin)
         if nargin > 2
            dfrmt = varargin{1};
         else
            dfrmt = obj.idataformat;
         end
         
         if ~isempty(dfrmt)
            d = impqlpermute(d, dfrmt, obj.irawformat);
         end 
      end
         
      
      function i = dataformatpos(obj, frmt)
         %
         % i = dataformatpos(obj, frmt)
         % 
         % description: 
         %    returns position of the format names in frmt in obj.idataformat
         
         [~, i] = ismember(frmt, obj.idataformat);
         %if length(frmt) ~= length(p)
         %   warning('ImageInfo: some image dimensions in %s are not found in %s', frmt, obj.iformat);
         %end 
      end
      
      function i = rawformatpos(obj, frmt)
         [~, i] = ismember(frmt, obj.irawformat);
      end
      
       function i = cellformatpos(obj, frmt)
         [~, i] = ismember(frmt, obj.icellformat);
       end

       

       function obj = renameDataFormat(obj, oldlab, newlab) 
          %
          % obj = renameDataFormat(obj, oldlab, newlab) 
          %
          % description:
          %    exchanges names of the dimensions oldlab(i) -> newlab(i)
          %    
          % input:
          %    oldlab     the format dimensions to change
          %    newlab     the new format dimensions 
           
          % check for conflict
          if any(ismember(setdiff(obj.idataformat, oldlab), newlab)) || length(oldlab) ~= length(newlab)
             error('ImageInfo: renameFormat: inconsistent reformatting of format %s from %s to %s', obj.idataformat, oldlab, newlab)
          end
                    
          ids = ismember(obj.idataformat, oldlab);
          obj.idataformat(ids) = newlab;
          obj.pqlctsizeFromFormatAndSize();
       end
       
       
       function obj = setDataFormat(obj, newformat)
          %
          % obj = setDataFormat(obj, newformat) 
          %
          % description:
          %    sets the data format to newformat and checks for consistency
          %    
          % input:
          %    newformat  the new data format
  
          % check for consistency
          if length(obj.idataformat) ~= length(newformat) || ~isempty(setdiff(obj.idataformat, newformat))
             error('ImageInfo: setFormat: inconsistent change of format fro %s to %s', obj.idataformat, newformat)
          end

          if isempty(obj.irawformat)
             obj.irawformat = obj.idataformat;
          end
          
          obj.idatasize   = obj.datasize(newformat);
          obj.idataformat = newformat;
       end

       
       function ids = seriesAsignIds(obj, varargin)
          %
          % ids = seriesAsignIds(obj, seriesid)
          %
          % decription:
          %     returns cell array such that data(ids{:}) can be assinged the series 
          %     data given the format specifications
          %
          % input:
          %     seriesid    (optional) the series number (=obj.iseries)
          %
          % output:
          %     ids         cell such that data(ids{:}) can be assinged the series data
          
          if nargin > 1
             s = varargin{1};
          else
             s = obj.iseries;
          end
          
          sd = obj.iseriesformat;
 
          ids = repmat({':'}, 1, obj.datadim);
          if isempty(sd)
             return
          end
          i = find(obj.idataformat == sd, 1, 'first');
          ids{i} = s;
       end

   end
end