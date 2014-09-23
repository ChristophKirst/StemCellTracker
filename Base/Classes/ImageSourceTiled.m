classdef ImageSourceTiled < ImageSource
   %
   % ImageSourceTiled class represents tiled image data
   %
   % description: organizes mapping from images to tiles in space
   % 

   properties 
      isource      = [];               % image source
      
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
         
         obj.isource = source;
      end
         
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% routines required by alignment 
      
      function id = tile2cellid(obj, id)
         per = imuvwformat2permute(obj.itileformat, 'uvw');
         id = imindpermute(obj.tilesize, per, id);
      end

      function img = getTile(obj, id)
         if obj.icachetiles
            imgs = obj.getTiles();
            img = imgs{obj.tile2cellid(id)};
         else
            % map id to id of the tileformatted data
            img = obj.isource.data(obj.tile2cellid(id));
         end
      end

      function imgs = getTiles(obj, varargin)
         if obj.icachetiles && ~isempty(obj.itiles)
            imgs = obj.itiles;
            if nargin > 1
               imgs =imgs(varargin{1});
            end
         else
            if nargin > 1
               ids = varargin{1};
               nids = numel(ids);
               imgs = cell(1,nids);
               for i = 1:nids
                  imgs{i} = obj.getTile(ids(i));
               end
               imgs = reshape(imgs, size(ids));
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
      
      function imgs = tiles(obj, varargin)
         imgs = obj.getTiles(varargin{:});
      end
      
      function ts = tilesizes(obj)
         ts = obj.getTileSizes();
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
      
      function n = ntiles(obj)
         n = prod(obj.tilesize);
      end
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% data access
      
      function d = getData(obj, varargin)
         if nargin > 0 && isnumeric(varargin{1})
            d = obj.isource.data(obj.tile2cellid(id));
         else
            d = obj.isource.data(varargin{:});
         end
      end
   
      function si = datasize(obj, varargin)
         si = obj.iinfo.datasize;
      end
      
      function ci = cellsize(~) 
         ci = 1;
      end


      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% utility
      
      function bkg = backgroundFromMinOfTiles(obj,varargin)
         bkg = backgroundFromMin(obj, varargin{:});
      end
      
      
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % info / visulaization
  
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
 

      function istr = infoString(obj)
         %istr = infoString@ImageSource(obj, 'Tiled');
         istr = 'ImageSource: Tiled';
         istr = [istr, '\ntileformat:     ', var2char(obj.itileformat)];
         istr = [istr, '\ntileshape:      ', var2char(obj.itileshape)];
         istr = [istr, '\ncachetiles:     ', var2char(obj.icachetiles)];  
         istr = [istr, '\n', obj.isource.infoString];    
      end
      
   end
      
   
end