classdef ImageSourceTiled < ImageSource
   %
   % ImageSourceTiling class represents a tiled image
   % 

   properties 
      isource      = [];               % image source
      ialignment   = ImageAlignment;   % image alignment
      
      itileformat = '';                % format of the tiling, usefull to permute tiling in correct way, to get tiles imuvwpermute(reshape(tiles, tileshape), tileformat, 'uvw') is used
      itileshape  = [];                % shape of the tiling  
      
      icachetiles = 1;                 % tile caching;
      itiles      = [];                % tile cache
   end

   methods
      function obj = ImageSourceTiled(varargin) % constructor
         %
         % ImageSourceTiled()
         % ImageSourceTiled(imagesourcetiling)
         % ImageSourceTiled(...,fieldname, fieldvalue,...)
         %

         if nargin == 0
            return
         elseif nargin >= 1
            if isa(varargin{1}, 'ImageSourceTiled') %% copy constructor
               obj = copy(varargin{1});
            elseif isa(varargin{1}, 'ImageSource')
               obj = obj.fromImageSource(varargin{:});
            %else
               %error('%s: invalid constructor input, expects char at position %g',class(obj), 1);
            %end
            %elseif nargin == 2 && isa(varargin{1}, 'ImageSource')  %% 
            %   obj.fromImageSource(varargin{:});
            else
               for i = 1:2:nargin % constructor from arguments
                  if ~ischar(varargin{i})
                     error('%s: invalid constructor input, expects char at position %g',class(obj), i);
                  end
                  if isprop(obj, lower(varargin{i}))
                     obj.(lower(varargin{i})) = varargin{i+1};
                  else
                     warning('%s: unknown property name: %s ', class(obj), varargin{i})
                  end
               end
            end
         end
      end
      
      
      function obj = fromImageSource(obj, source, varargin)
         param = parseParameter(varargin);
         ts = getParameter(param, 'tileshape', source.cellsize);
         tf = 'uvw';
         tf = getParameter(param, 'tileformat', tf(1:length(ts)));

         obj.itileshape  = ts;
         obj.itileformat = tf;
         
         obj.ialignment = ImageAlignment(ts);
         obj.isource = source;
      end
         
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% routines required by alignment 
      
      function id = tile2cellid(obj, id)
         per = imuvwformat2permute(obj.itileformat, 'uvw');
         id = imindpermute(obj.tilesize, per, id);
      end

      function img = getTile(obj, id)
         % map id to id of the tileformatted data   
         img = obj.isource.data(obj.tile2cellid(id));
      end

      function imgs = getTiles(obj)
         if obj.icachetiles && ~isempty(obj.itiles)
            imgs = obj.itiles;
         else
            imgs = obj.isource.celldata;
            
            if ~isempty(obj.itileformat)
               imgs = reshape(imgs, obj.itileshape);
               imgs = imuvwpermute(imgs, obj.itileformat, 'uvw');
            end
            
            if obj.icachetiles
               obj.itiles = imgs;
            end
         end
      end
      
      function obj = clearCache(obj)
         clearCache@ImageSource(obj);
         obj.itiles = [];
      end
      
      function si = getTileSizes(obj)
         ti = obj.isource.datasize;
         
         if ~isempty(obj.itileformat)
            ci = obj.tilesize;
            ci = ci(imuvwformat2permute(obj.itileformat, 'uvw'));
         else
            ci = obj.isource.cellsize;
         end
         
         if length(ci) == 1
            ci(2) = 1;
         end
         si = repmat({ti}, ci);
      end
      
      
      function setTileFormat(obj, tfrmt)
         obj.itileformat = tfrmt;
      end
      
      
      function tf = tileformat(obj)
         tf = obj.itileformat;
      end
      
      %%% same as above 
      function img = tile(obj, id)
         img = obj.getTile(id);
      end
      
      function imgs = tiles(obj)
         imgs = obj.getTiles();
      end
      
      function imgs = tileSizes(obj)
         imgs = obj.getTileSizes();
      end
      
      function tsi = tilesize(obj)
         if ~isempty(obj.itileshape)
            tsi = obj.itileshape;
            
            if ~isempty(obj.itileformat)
               tsi = tsi(imuvwformat2permute(obj.itileformat, 'uvw'));
            end
         else
            tsi = obj.isource.cellsize;
         end
      end
      
      function td = tiledim(obj)
         td = length(obj.tilesize);
      end
      

      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% data access

      function d = getData(obj, varargin)
         d = obj.ialignment.stitch(obj, varargin{:}); %if cache is on this will be cached automatically after first run when using the data routine
      end
   

%       function d = subAlignment(obj)
%          
%          
%       end


      function si = datasize(obj)
         si = obj.ialignment.iasize;
      end
      
      function ci = cellsize(~) 
         ci = 1;
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
      
 
      function plotAlignedImages(obj)
         %
         % plotAlignedImages(obj)
         %
         % description:
         %    visualizes the alignment using plotAlignedImages
         %
         % See also: plotAlignedImages
                  
         plotAlignedImages(obj.getTiles, obj.imageShifts, 'colors', obj.tiledim);
      end
 
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % infoString
      function istr = infoString(obj)
         istr = infoString@ImageSource(obj, 'File');
         istr = [istr, '\nfilename:   ', obj.ifilename];
         istr = [istr, '\nreader:     ', func2str(obj.ireader)];
         istr = [istr, '\ninforeader: ', func2str(obj.iinforeader)]; 
      end
      
   end
      
   
end