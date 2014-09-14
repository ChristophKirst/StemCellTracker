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
     
      irawsize   = [];               % size of the raw data, [] = isize
      irawformat = '';               % format of the raw image data, when read from a file, '' = iformat
  
      idataclass = '';               % class of the image
            
      icolor = {};                   % colors for the channels
      
      iseries = 1;                   % (optional) series number or ids if image if from a larger file 
      iseriesformat = '';            % (optional) dimension the series is linked to 

      inimages = [];                 % (optional) number of total 2d images as counted via loci tools

      imetadata = [];                % (optional) meta data
      
      iscale = [];                   % (optional) spatial scale of image in pixel per spatial unit
      iunit = '';                    % (optional) spatial unit
      
      iname  = [];                   % (optional) name of the image 
      
      %ikeys = ImageKeys;           % translation between shortcut key (e.g. 'dapi') and a subset of the data (e.g. 'c' -> 1)
   end
   
   methods
            
      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% init
      

      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% sizes
      
      function d = datadim(obj)
         d = length(obj.idataformat);
      end
      
      function d = celldim(obj)
         d = length(obj.icellformat);
      end
      
      
      
      function obj = setDataSize(obj, newsize) %using this routine is usally not a good idea
         obj.idatasize = newsize;
         obj.pqlctsizeFromFormatAndSize();
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


      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% formatting
      
      function obj = pqlctsizeFromFormatAndSize(obj)
         frmt = 'pqlct';  
         pos = arrayfun(@(x) find(x == frmt,1), obj.idataformat);
         obj.idatasizePQLCT = ones(1, 5);
         obj.idatasizePQLCT(pos) = obj.idatasize;
      end
      
      function d = raw2dataformat(obj, d, varargin)
         if nargin > 2
            rfrmt = varargin{1};
         else
            rfrmt = obj.irawformat;
         end
         
         if ~isempty(rfrmt)
            d = impqlpermute(d, rfrmt, obj.idataformat);
         end
      end
      
      function s = datasize(obj, varargin)
         if nargin == 1
            s = obj.idatasize;
         else
            s = obj.idatasize(obj.dataformatpos(varargin{1}));
         end
      end
      
      
      
      
      function i = dataformatpos(obj, frmt)
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

       
       function si = dataformatsize(obj, frmt)
          i = obj.dataformatpos(frmt);
          si = obj.idatasize;
          si = si(i);
       end
       
       function si = rawformatsize(obj, frmt)
          i = obj.rawformatpos(frmt);
          si = obj.irawsize;
          si = si(i);
       end
       
       
       % exchanges names of the dimensions
       function obj = renameDataFormat(obj, oldlab, newlab) 
          % check for conflict
          if any(ismember(setdiff(obj.idataformat, oldlab), newlab)) || length(oldlab) ~= length(newlab)
             error('ImageInfo: renameFormat: inconsistent reformatting of format %s from %s to %s', obj.idataformat, oldlab, newlab)
          end
                    
          ids = ismember(obj.idataformat, oldlab);
          obj.idataformat(ids) = newlab;
          obj.pqlctsizeFromFormatAndSize();
       end
       
       
       function obj = setDataFormat(obj, newformat)
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

       
      
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%% data asignments / permutations
       
       function ids = seriesAsignIds(obj, varargin)
          % returns cell array such that data(ids{:}) can be assinged the series data
          
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

       function dat = permuteData(obj, dat, varargin)
          if nargin > 2
             frmt = varargin{1};
          else
             frmt = obj.idataformat;
          end

          dat = impqlpermute(dat, obj.irawformat, frmt);
       end

   end
end