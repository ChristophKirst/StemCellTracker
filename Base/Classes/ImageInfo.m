classdef ImageInfo < matlab.mixin.Copyable
   
   % class that manages information about a single image in a file
   % for multiple series in a file choose ImageSourceList
   
   % notes:
   %    the info class is based on the info needed to read images via imread_bf / loci tools
   %    data read from files always linked to the pqlct dims to simplify array manipulations / reading etc
   %    series can be linked to cell or image dims
   %    tags can be linked to cell or image dims
   
   properties
      isize      = [];             % size of the image as returned by .data routine 
      iformat    = '';             % format of the image as returned by .data routine
      
      isizePQLCT = [0,0,0,0,0];    % size of the image for full pqlct dimensions

      icellsize   = 1;             % size of the cell structure as returend when using .celldata routine
      icellformat = '';            % format of the cell structure (usually fixed subset of uvwrs)
      
      %icellsizeUVWRS = [0,0,0,0,0] % size of the cell for full uvwrs dimensions 
      
      irawformat = '';             % format of the raw image data when read from a file, '' = iformat
      %irawcellformat = '';        % format of the raw cell structure , '' = icellformat

      iclass = '';                 % class of the image
      
      icolor = {};                 % colors for the channels
      
      iseries = 1;                 % (optional) series number or ids if image if from a larger file 
      iseriesdim = '';             % (optional) dimension the series is linked to 

      inimages = [];               % (optional) number of total 2d images as counted via loci tools

      imetadata = [];              % (optional) meta data
      
      iscale = [];                 % (optional) spatial scale of image in pixel per spatial unit
      iunit = '';                  % (optional) spatial unit
      
      iname  = [];                 % (optional) name of the image 
      
      %ikeys = ImageKeys;           % translation between shortcut key (e.g. 'dapi') and a subset of the data (e.g. 'c' -> 1)
   end
   
   methods
            
      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% init
      
      function obj = initFromBf(obj)
         obj.irawformat = obj.iformat;
         obj.icellsize = length(obj.iseries);
      end
      

      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% sizes
      
      function d = dim(obj)
         d = length(obj.iformat);
      end
      
      function d = celldim(obj)
         d = length(obj.icellformat);
      end
      
      function obj = setSize(obj, newsize) %using this routine is usally not a good idea
         obj.isize = newsize;
         obj.pqlctsizeFromFormatAndSize();
      end
      
      
 
      function p = sizeP(obj)
         p = obj.isizePQLCT(1);
      end
      function q = sizeQ(obj)
         q = obj.isizePQLCT(2);
      end
      function l = sizeL(obj)
         l = obj.isizePQLCT(3);
      end
      function c = sizeC(obj)
         c = obj.isizePQLCT(4);
      end
      function t = sizeT(obj)
         t = obj.isizePQLCT(5);
      end
      
%       function p = sizeU(obj)
%          p = obj.isizeUVWRS(1);
%       end
%       function q = sizeV(obj)
%          q = obj.isizeUVWRS(2);
%       end
%       function l = sizeW(obj)
%          l = obj.isizeUVWRS(3);
%       end
%       function c = sizeR(obj)
%          c = obj.isizeUVWRS(4);
%       end
%       function t = sizeS(obj)
%          t = obj.isizeUVWRS(5);
%       end


      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% formatting
      
      function obj = pqlctsizeFromFormatAndSize(obj)
         frmt = 'pqlct';  
         pos = arrayfun(@(x) find(x == frmt,1), obj.iformat);
         obj.isizePQLCT = ones(1, 5);
         obj.isizePQLCT(pos) = obj.isize;
      end
      
      function d = rawdata2data(obj, d)
         d = impqlpermute(d, obj.idataformat, obj.iformat);
      end
      
      function s = size(obj, varargin)
         if nargin == 1
            s = obj.isize;
         else
            s = obj.isize(obj.formatpos(varargin{1}));
         end
      end
      
      function i = formatpos(obj, frmt)
         [~, i] = ismember(frmt, obj.iformat);
         %if length(frmt) ~= length(p)
         %   warning('ImageInfo: some image dimensions in %s are not found in %s', frmt, obj.iformat);
         %end 
      end
      
      function i = dataformatpos(obj, frmt)
         [~, i] = ismember(frmt, obj.idataformat);
      end
      
       function i = cellformatpos(obj, frmt)
         [~, i] = ismember(frmt, obj.icellformat);
       end
         
       
       % exchanges names of the dimensions
       function obj = renameFormat(obj, oldlab, newlab) 
          % check for conflict
          if any(ismember(setdiff(obj.iformat, oldlab), newlab)) || length(oldlab) ~= length(newlab)
             error('ImageInfo: renameFormat: inconsistent reformatting of format %s from %s to %s', obj.iformat, oldlab, newlab)
          end
                    
          ids = ismember(obj.iformat, oldlab);
          obj.iformat(ids) = newlab;
          obj.pqlctsizeFromFormatAndSize();
       end
       
       
       function obj = permuteFormat(obj, newformat)
          % check for consistency
          if length(obj.iformat) ~= length(newformat) || ~isempty(setdiff(obj.iformat, newformat))
             error('ImageInfo: setFormat: inconsistent change of format fro %s to %s', obj.iformat, newformat)
          end

          if isempty(obj.idataformat)
             obj.idataformat = obj.iformat;
          end
          
          obj.isize   = obj.size(newformat);
          obj.iformat = newformat;
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
          
          sd = obj.iseriesdim;
 
          ids = repmat({':'}, 1, obj.dim);
          if isempty(sd)
             return
          end
          i = find(obj.iformat == sd, 1, 'first');
          ids{i} = s;
       end

       function dat = permuteData(obj, dat, varargin)
          if nargin > 2
             frmt = varargin{1};
          else
             frmt = obj.iformat;
          end
          
          dat = impqlpermute(dat, obj.idataformat, frmt);
       end

       function img = permuteImage(obj, img, frmt)
          img = impqlpermute(img, obj.iformat, frmt);
       end
      
   end
end