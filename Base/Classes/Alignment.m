classdef Alignment < ImageSource
   %
   % Alignment class representing alignment/tiling data
   %

   properties
      apairs = [];       % pairs of overlapping images (array of AlignmentPair classes)
      anodes = [];       % ids of images used for this alignment
            
      asource = [];      % image source. assumed to have routines cell(id), cellSize(id) and data(id)
      aformat = [];      % coordinates in the cell dims of the source that contribute to the tiling

      aposition = [1,1]; % position of the first node (i.e. if connected component of a larger tiling)
      ashifts = {};      % image shifts assuming first node at (1,1,...)

      aerror  = Inf;     % alignment error
   end
  
   methods

      function obj = Alignment(varargin)
         %
         % Alignment()
         % Alignment(cellarray)
         % Alignment(struct)
         % Alignment(..., fieldname, fieldvalue, ...)
         %
         if nargin >= 1 
            if isa(varargin{1}, 'Alignment') %% copy constructor
               obj = copy(varargin{1});
            elseif isa(varargin{1}, 'AlignmentPair')
               obj.fromAlignmentPair(varargin{1});
            elseif isa(varargin{1}, 'cell')
               obj.fromCell(varargin{1});
            elseif isa(varargin{1}, 'struct')
               obj.fromStruct(varargin{1});
            elseif isa(varargin{1}, 'ImageSource')
               obj.fromImageSource(varargin{:});  
            else
               error('%s: not valid arguments for constructor', class(obj));
            end
         else
            obj.fromParameter(varargin);
         end
         
         obj.setCaching(false);  % caching the stitched image is usually not a good idea if large
      end

      function obj = fromParameter(obj, varargin)
         obj = classFromParameter(obj, 'a', varargin);
      end

      function obj = fromCell(obj, ca)
         %
         % obj = fromCell(obj, ca)
         %
         % descritpion:
         %   constructs Alignment class from prealigned cell array ca
         %
    
         si = size(ca);
         si = padright(si, 3, 1);

         np = si(1) * si(2) * (si(3)-1) + si(1) * (si(2)-1) * si(3) + (si(1)-1) * si(2) * si(3);
         p(np) = AlignmentPair();

         i = 1;
         for x = 1:si(1)
            for y = 1:si(2)
               for z = 1:si(3)
                  if x < si(1)
                     p(i).from        = ca(x,y,z);
                     p(i).to          = ca(x+1,y,z);
                     p(i).orientation = 1;
                     i = i + 1;
                  end
                  if y < si(2)
                     p(i).from = ca(x,y,z);
                     p(i).to   = ca(x,y+1,z);
                     p(i).orientation = 2;
                     i = i + 1;
                  end
                  if z < si(3)
                     p(i).from = ca(x,y,z);
                     p(i).to   = ca(x,y,z+1);
                     p(i).orientation = 3;
                     i = i + 1;
                  end
               end
            end
         end
         
         obj.apairs = p;
         obj.anodes = [ca{:}];
      end

      function obj = fromCellSize(obj, si)
         %
         % obj = fromCellSize(obj, ca)
         %
         % descritpion:
         %   constructs Alignment class from prealigned cell array ca
         %
    
         %si = size(ca);
         si = padright(si, 3, 1);

         np = si(1) * si(2) * (si(3)-1) + si(1) * (si(2)-1) * si(3) + (si(1)-1) * si(2) * si(3);
         p(np) = AlignmentPair();

         i = 1;
         for x = 1:si(1)
            for y = 1:si(2)
               for z = 1:si(3)
                  if x < si(1)
                     p(i).from        = imsub2ind(si, [x,y,z]);
                     p(i).to          = imsub2ind(si, [x+1,y,z]);
                     p(i).orientation = 1;
                     i = i + 1;
                  end
                  if y < si(2)
                     p(i).from = imsub2ind(si, [x,y,z]);
                     p(i).to   = imsub2ind(si, [x,y+1,z]);
                     p(i).orientation = 2;
                     i = i + 1;
                  end
                  if z < si(3)
                     p(i).from = imsub2ind(si, [x,y,z]);
                     p(i).to = imsub2ind(si, [x,y,z+1]);
                     p(i).orientation = 3;
                     i = i + 1;
                  end
               end
            end
         end
         
         obj.apairs = p;
         obj.anodes = 1:prod(si);
      end

      function obj = fromStruct(obj, st)
         %
         % obj = fromStruct(obj, st)
         %
         % descritpion:
         %   sets pair info from a pair of images given by a prealigned cell array 
         %
         % input:
         %   st  struct with .from, .to and additional information such as shift, orientation, quality
         %
         
         obj.apairs = AlignmentPair(st);
         obj.anodes = obj.nodeIndexFromPairs();
      end
      
      function obj = fromAlignmentPair(obj, ap)
         %
         % obj = fromAlignmentPair(obj, ap)
         %
         % descritpion:
         %   sets pair info from an array of AlignmentPair classes 
         %
         % input:
         %   ap  array of AlignmentPair classes 

         obj.apairs = AlignmentPair(ap);  % deep copy
         obj.anodes = obj.nodeIndexFromPairs();
      end
      
      function obj = fromImageSource(obj, is, frmt)
         %
         % obj = fromImageSource(obj, is)
         %
         % descritpion:
         %   intitializes class from an ImageSource class
         %
         % input:
         %   is    ImageSource class
         %   fmrt  (optional) the coordinates for the tiles
         
         if nargin < 3
            obj.aformat = is.cellFormat;
         else
            pos = imfrmtPosition(is.cellFormat, frmt);
            frmt = is.cellFormat;             
            obj.aformat = frmt(pos);
         end    
         % check if non-tile dims are singletons
         
         noFrmt = imfrmtRemoveFormat(is.cellFormat, obj.aformat);
         if any(is.cellSize(noFrmt) > 1)
            error('%s: fromImageSource: not all non tiling dimensions are singletons in source!');
         end
         
         obj.asource = is;
         obj.fromCellSize(is.cellSize(obj.aformat)); 
         
         obj.ipreviewscale = is.previewScale;
      end
      
      
      function  obj = fromImageSourceBF(obj, is, varargin)
         
         % find grid of nodes form the positions assuming a grid like structure
         
         param = parseParameter(varargin);
         
         pos = getParameter(param, 'positions', []);
         if isempty(pos)
            pos = is.stagePositions;
         end

         isiz = [is.ireader.getSizeX; is.ireader.getSizeY];
         vsiz = imreadBFVoxelSize(is.ireader, 'S', 1, varargin);
         vsiz = vsiz{1};

         isiz = isiz .* vsiz(1:2); 

         ni = length(pos);
         pos = [pos{:}];
         posx = pos(1,:)'; posy = pos(2,:)';

         % x clustering
         
         [clx, idx] = clusterMeanShift(posy, isiz(1)/2)
         ntx = length(clx)

         if mod(ni, ntx) ~= 0
            error('%s: fromImageSourceAndStagePositions: cannot infer gird!', class(obj));
         end
         
         % y clustering
         
         [cly, idy] = clusterMeanShift(posx, isiz(2)/2);
         nty = length(cly);

%          if mod(ni, nty) ~= 0
%             error('%s: fromImageSourceAndStagePositions: cannot infer gird!', class(obj));
%          end
         
         % z clustering
         if size(pos,1) > 2
            posz = pos(3,:)';
            
            isz = is.ireader.getSizeZ * vsiz(3);
            [clz, idz] = clusterMeanShift(posz, isz/2);

            ntz = length(clz);

%             if mod(ni, ntz) ~= 0
%                error('%s: fromImageSourceAndStagePositions: cannot infer gird!', class(obj));
%             end
         else
            ntz = 1;
            idz = ones(1,ni);
         end
         
%          if ntz*nty*ntz ~= ni
%             error('%s: fromImageSourceAndStagePositions: cannot infer gird, dimensions mistmatch!', class(obj));
%          end

         % fill grid
         gr = cell(ntx, nty, ntz);
         ids = sub2ind([ntx, nty, ntz], idx, idy, idz);
         gr(ids) = 1:ni;
         
         % 
         obj = obj.fromCell(gr);
         
         obj.ipreviewscale = is.previewScale;
      end

      
      function  obj = fromStagePositionsAndVoxelSize(obj, pos, vsiz)
         %
         % obj = fromStagePositions(obj, pos, vs)
         %
         % initializes shifts from positions
         
         % intialize shifts directly form stage positions
         if nargin == 2 && isa(pos, 'ImagesouceBF')
            obj.asource = pos;
 
            vsiz = imreadBFVoxelSize(pos.ireader, 'S', 1, varargin);
            pos = pos.stagePositions;
            
            obj.anodes = 1:length(pos);
         end

         % to pixel coords
         
         vsiz = vsiz(1:length(pos{1}));
         pos = cellfunc(@(x) round(x ./ vsiz), pos);
         
         pos = cellfunc(@(x) x(:)', pos);
         obj.ashifts = cellfunc(@(x) x - pos{1}, pos);

         obj.initializeFromShifts;
      end
      

      function obj = initializeFromShifts(obj)
         % intializes the data size and format info using the info in pairs
         
         obj.irawdatasize = obj.absoluteSize;
         frmt = 'XYZ';
         frmt = frmt(1:length(obj.irawdatasize));
         obj.irawdataformat = frmt;
         
         obj.irawcellsize = 1;
         obj.irawcellformat = 'S';
         
         obj.initializeCellDataSizeAndFormatFromRaw;
      end
         

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % basics 
      
      function d = dim(obj)
         %
         % d = dim(obj)
         %
         % descritpion:
         %   number of dimensions of the alignment shifts

         if ~isempty(obj.apairs)
            d = obj.apairs(1).dim();
         else
            d = 2;
         end
      end
      
      function n = nPairs(obj)
         %
         % n = npairs(obj)
         %
         n = cellfun(@length, {obj.apairs});
      end
      
      function n = nNodes(obj)
         %
         % n = nnodes(obj)
         %
         
         n = cellfun(@length, {obj.anodes});
      end
      
      function n = nodes(obj)
         n = obj.anodes;
      end
      
      function p = pairs(obj)
         p = obj.apairs;
      end
      
      function s = source(obj)
         s = obj.asource;
      end
      
      
      function idx = cellIndexFromNodeIndex(obj, id)
         idx = obj.anodes(id);
      end
      
      function idx = nodeIndexFromCellIndex(obj, varargin)
         id = obj.asource.cellIndex(varargin{:});
         [~, idx] = ismember(id, obj.anodes);
      end
      
      function id = nodeIndexFromPairs(obj)
         %
         % id = indexFromPairs(obj)
         %
         % descritpion:
         %   returns all ids used in pairs
         
         id = unique([[obj.pairs.from], [obj.pairs.to]]);
      end
 
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% underlying image source
      
      function ts = sourceCellSize(obj, varargin)
         ts = obj.asource.cellSize(varargin{:});
      end
      
      function tf = sourceCellFormat(obj)
         tf = obj.asource.cellFormat;
      end

      function cd = sourceCell(obj, varargin)
         cd = obj.asource.cell(varargin{:});
      end
      
      function d = sourceData(obj, varargin)
         d = obj.asource.data(varargin{:});
      end

      function is = sourceDataSize(obj, varargin)
         is = obj(1).asource.dataSize(varargin{:});
      end
      
      function is = sourceDataFormat(obj, varargin)
         is = obj.asource.dataFormat(varargin{:});
      end
      
      function d = nodeData(obj, varargin)
         ids = obj.cellIndexFromNodeIndex(varargin{:});
         d = obj.asource.data(ids);
      end
         
      function d = nodeCell(obj, varargin)
         ids = obj.cellIndexFromNodeIndex(varargin{:});
         d = obj.asource.cell(ids);
      end
      
      function obj = sourceClearCache(obj)
         obj = obj.asource.clearCache;
      end
      

  
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % images / tiles / positions

      function p = position(obj)
         p = obj.aposition;
      end
      
      function obj = setPosition(obj, o)
         obj.aposition = o;
      end
       
      function p = absolutePosition(obj)
         ap = absoluteShiftsAndSize(obj.imagePositions, obj.imageSizes);
         k = 1;
         n = length(obj);
         p = cell(1, n);
         for i = 1:n
            p{i} = ap{k} + 1;
            k = k + obj(i).nNodes;
         end
      end
      
      function si = absoluteSize(obj)
         [~,si] = absoluteShiftsAndSize(obj.imagePositions, obj.imageSizes);
      end
      
      
      function imgs = images(obj, varargin)
         nodes = [obj.anodes];
         if nargin == 2
            nodes = nodes(varargin{:});
         end   
         source = obj.asource; 
         imgs = source.cell(nodes);
         imgs = imgs(:);
      end
      
      
      function imgs = imagesResmaple(obj, scale, varargin)
         nodes = [obj.anodes];
         if nargin == 2
            nodes = nodes(varargin{:});
         end   
         source = obj.asource; 
         imgs = source.cellResample(scale, nodes);
         imgs = imgs(:);
      end

      function isizes = imageSizes(obj, varargin)
         nodes = [obj.anodes];
         if nargin == 2
            nodes = nodes(varargin);
         end 
         isizes = repmat({obj.sourceDataSize}, 1, length(nodes)); 
      end

      
      function shifts = imageShifts(obj, varargin)
         %
         % shifts = imageShifts()
         %
         % description:
         %    returns shifts for the images w.r.t to node 1, assumes ashifts are set
         %
         % See also: align

         shifts  = vertcat(obj.ashifts);
         if isempty(shifts)
            error('%s: imageShifts: shifts not set, align images first !', class(obj));
         end

         shifts = shifts(varargin{:})';
      end
      

      function [pos, si] = absoluteImageShifts(obj)
         %
         % obj = absoluteImageShifts(obj, varargin)
         %
         % description:
         %     returns shifts so that all image coords are positive

         [pos, si] = absoluteShiftsAndSize(obj.imagePositions, obj.imageSizes);
      end
      

      function ipos = imagePositions(obj, varargin)
         %
         % ipos = imagePositions()
         %
         % description:
         %    returns 'physical' positions of the images given the shifts and position
         %
         % See also: align

         ipos = {obj.ashifts};

         if any(cellfun(@isempty, ipos))
            error('%s: imageShifts: shifts not set, align images first !', class(obj));
         end
         
         for i = 1:length(ipos)
            o = obj(i).aposition;
            ipos{i} = cellfunc(@(x) x + o, ipos{i});
         end
         
         ipos = vertcat(ipos{:});
         
         ipos = ipos(varargin{:});
      end
            
      function pos = absoluteImagePositions(obj)
         %
         % obj = absoluteImagePositions(obj, varargin)
         %
         % description:
         %     returns positions so that all image coords are positive

         pos = absoluteShiftsAndSize(obj.imagePositions, obj.imageSizes);
         pos = cellfunc(@(x) x + 1, pos);
      end
      

      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% data access
      
       % raw data access routine to be overwritten by super class 
      function d = getRawData(obj, range)
         %
         % d = getRawData(obj, varargin)
         %
         % description:
         %     obtain raw image data in the format given by obj.rawformat
         %     data specs use raw formats and sizes

         cid = obj.rawCellIndexFromRawVarargin(range);
         if length(cid) > 1
            error('%s: getRawData: cell dimensinos are not specified to be singeltons!', class(obj));
         end
         
         d = obj.stitch;
         d = imfrmtDataSubset(d, obj.rawDataFormat, range);
      end

      % raw cell data access routine to be overwritten by super class 
      function d = getRawCellData(obj, range)
         %
         % d = getRawData(obj, varargin)
         %
         % description:
         %     obtain raw cell data in the format given by obj.rawformat
         %     data specs use raw formats and sizes
                  
         d = {obj.getRawData(range)}; % only single stitched image in here
      end
      
      
      function [d, sh] = dataExtract(obj, roi, varargin)
         %
         % d = dataExtract(obj, roi)
         %
         % description:
         %    extract a subset of the data given the spatial roi 
         %
         
         if obj.irawcache
            if isempty(obj.irawcelldata)
               obj.stitch(varargin{:});
            end
            
            [d, sh] = roi.extractData(obj.data);
            
         else % serial processing

            % we dont want to stitch a huge image for a small roi 
            % -> find relevant ids and pairs first
            [ids, shids] = obj.nodesFromROI(roi);

            ash = obj.imageShifts;
            ash = ash(shids);
            tsi = obj.imageSizes;
            tsi = tsi(shids);

            st = stitchImages(obj.sourceCell(ids), ash, varargin{:});
            
            as = absoluteShiftsAndSize(ash, tsi);
            
            % correct roi
            r = roi.copy;
            r.shift(as{1}-ash{1});

            % extract
            [d, sh] = r.extractData(st);
            
            sh = sh + ash{1}-as{1};
         end
      end
      
      

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % alignment
      
      %%% routines modifying the class
      
      function obj = removeLowQualityPairs(obj, thq)
         %
         % obj = removeLowQualityPairs(obj, thq)
         %
         % descritpion:
         %   removes pairs whose quality is less than thq
         %
         % input:
         %   thq   quality threshold 
         
         q = [obj.apairs.quality];
         obj.apairs = obj.apairs(q >= thq); 
      end
      
      function obj = reducePairs(obj)
         %
         % obj = reduce(obj)
         %
         % description: 
         %     removes all paris not connecting two nodes in obj.nodes

         iids = obj.anodes;
          
         id  = [obj.apairs.from];
         ids = arrayfun(@(x) any(x==iids), id);
                  
         id  = [obj.apairs.to];
         ids = and(ids, arrayfun(@(x) any(x==iids), id));
         
         obj.apairs = obj.apairs(ids);
      end
      
      function obj = removeBackgroundNodes(obj, thb)
         %
         % comp = removeBackgroundNodes(obj, thb)
         %
         % descritpion:
         %   given an array of alignments remove the ones in the node list
         %   for which there is no signal in the images
         %
         % input:
         %   thb    threshold for background intensity
    
         nds = obj.anodes;
         n = obj.nNodes;
         rmn = zeros(1,n);
         for i = 1:n
            img = obj.asource.data(nds(i));
            rmn(i) = max(img(:)) < thb;
         end
         
         if any(rmn)
            obj.nNodes = nds(~rmn);
            obj.reducePairs();
         end
      end
      
      function obj = calculateOverlapQuality(obj, varargin)
         %
         % obj = overlapQuality(obj) 
         %
         % description:
         %    calculates operlap quality of the images
         %
         % See also: overlapQuality, overlapStatisticsImagePair
         pairs = obj.apairs;
         source = obj.asource;
         
         parfor p = 1:obj.nPairs
            pp = pairs(p).copy();
            pp.from = source.data(pp.from); %#ok<PFBNS>
            pp.to   = source.data(pp.to); 
            pairs(p).quality = overlapQuality(pp, varargin{:}); %#ok<PFBNS>
         end
         
         obj.apairs = pairs;
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
         
         if isempty(obj.asource)
            error('alignment: alignParis: no image data to align images, set property asource');
         end
         
         pairs = obj.apairs;
         source = obj.asource;

         parfor i = 1:obj.npairs
            p = AlignmentPair(pairs(i)); % deep copy
            p.from = source.data(p.from); %#ok<PFBNS>
            p.to   = source.data(p.to);
            
            p = alignImagePair(p, varargin{:}); %#ok<PFBNS>
           
            % we dont want to copy the data but keep the indices
            pairs(i).shift   = p.shift;
            pairs(i).quality = p.quality;
            pairs(i).aerror  = p.aerror;
         end
         
         obj.apairs = pairs;
      end
     
      function obj = optimizePairwiseShifts(obj)
         %
         % obj = optimizePairwiseShifts(obj)
         %
         % description:
         %    globally optimizes pairwise shifts 
         %
         % See also: optimizePairwiseShifts
         
         [obj.apairs, ic] = optimizePairwiseShifts(obj.pairs);
         ic = round(ic);
   
         % transform consistent shifts to shifts 
         shifts = num2cell(ic',2);

         sh = shifts{1};
         for i = 1:numel(shifts)
            shifts{i} = shifts{i} - sh;
         end
         obj.ashifts = shifts;
         
         obj.initializeFromShifts();
      end
      
      
      %%%  align and stitch
      
      function obj = align(obj, varargin)
         %
         % obj = align(obj, varargin)
         %
         % description:
         %    aligns images and sets new shifts
         %
         % See also: alignImages
         
         [~, pairs] = alignImages(obj.asource, 'pairs', obj.apairs, 'nodes', obj.anodes,  varargin{:});
         obj.apairs = pairs;
         
         if obj.nNodes == 1
            obj.ashifts = {zeros(1,obj.dim)};
            obj.initializeFromShifts();
            return
         end
         
         [~, ic] = optimizePairwiseShifts(obj.apairs);
         ic = round(ic);
   
         % transform consistent shifts to shifts 
         shifts = num2cell(ic',2);
         %var2char(shifts);

         sh = shifts{1};
         for i = 1:numel(shifts)
            shifts{i} = shifts{i} - sh;
         end
         obj.ashifts = shifts;
         
         obj.initializeFromShifts();
      end

      
      function st = stitch(obj, varargin)
         %
         % st = stitch(obj, param)
         %
         % description
         %     stitches the images using the alignment information and parameter params
         
         st  = stitchImages(obj.images, obj.imageShifts, varargin{:});
         
         if obj.irawcache  % cach if required
            obj.fromData(st);
         end
      end
      
      function st = stitchArray(obj, varargin)
         % stitches an array of alignement classes using thier different origins
         st  = stitchImages(obj.images, obj.imagePositions, varargin{:});
      end
      
            
      function st = stitchArrayResample(obj, scale, varargin)
         % stitches an array of alignement classes using their different origins 
         % resampling the data
         pos = obj.imagePositions;
         pos = cellfunc(@(x)round(x .* scale), pos);
         st  = stitchImages(obj.imagesResample(scale), pos, varargin{:});
      end
      
      
      
      function aerror = alignmentError(obj)
          if obj.nPairs == 0
              aerror = 0;
          else
            aerror = sum([obj.apairs.aerror]) / obj.nPairs;
          end
      end

      
  
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % connected components / origins

      function oq = overlapQuality(obj)
         oq = [obj.apairs.quality];
      end

      function comp = connectedComponents(obj, varargin)
         %
         % comp = connectedComponents(obj, varargin)
         %
         % descritpion:
         %   given quality factors of pairs returns cell array of Alignemnts with connected images
         %
         % input:
         %   param   parameter as in connectedAlignments
         %
         % See also: connectedAlignments
         
         obj.clearPreview;
         comp = connectedAlignments(obj, varargin{:});   
      end
 
      function ov = overlapSizePrimary(obj)
         %
         % ov = meanOverlapSizePrimary(obj)
         %
         % description:
         %    calculates man overlap size in primary directinos form the pairs

         nobj = length(obj);
         dm = obj(1).dim;
         ds = obj(1).sourceDataSize;
         ov = cell(1,3);
         
         for i = 1:nobj;
            or = [obj(i).apairs.orientation];
            ovl = obj(i).apairs.overlapSizePrimary(ds);

            for d = 1:dm
               ov{d} = [ov{d}, ovl(or == d)];
            end
         end

         ov = ov(1:dm);
      end
      
      function ov = meanOverlapSizePrimary(obj)
         %
         % ov = meanOverlapSizePrimary(obj)
         %
         % description:
         %    calculates man overlap size in primary directinos form the pairs

         ov = obj.overlapSizePrimary;
         ov = cellfun(@mean, ov);
      end
      
     
      
      function obj = alignPositions(obj, varargin)
         %
         % obj = alignPositions(obj, varargin)
         %
         % description:
         %    aligns origins of an array of Alignment classes originating form a common grid
         %    e.g. as obtained by connectedComponents
         %
         % input
         %    parameter struct with entries
         %             overlap.mean    values to assume for mean overlap ([] = automatic)
         %             overlap.scale   scaling of the mena overlap to get obtain the final overlap ([] = 1)     
         %
         % See also: connectedComponents

         nobj = length(obj);

         
         % overlap
         param = parseParameter(varargin);
         
         ovlm = getParameter(param, 'overlap.mean', []);
         if isempty(ovlm)
            ovlm = obj.meanOverlapSizePrimary;
         end
         
         ovls = getParameter(param, 'overlap.scale', 1);
         ovls = padright(ovls, length(ovlm), ovls);
         
         ovl = ovlm .* ovls;

         % assume the same source for all 
         isource = obj(1).source;
         %pos = imfrmtPermuation(isource.cellFormat, obj(1).aformat);
         pos = imfrmtPosition(isource.cellFormat, obj(1).aformat);
         
         isiz = obj.sourceDataSize;

         for i = 1:nobj
            nodes = obj(i).nodes;
            sub = isource.cellIndexToSubIndex(nodes(1));
            sub = sub(:, pos);
            
            obj(i).aposition = round((sub - 1) .* (isiz - ovl)) + 1;
         end
      end
      

      
      
      function algn = merge(obj)
         %
         % algn = merge(obj)
         %
         % description:
         %   merges the sub alignment classes to a single one
         
         
         algn = obj(1).copy;
         algn.clearPreview;
         
         algn.apairs = [obj.apairs];
         algn.anodes = [obj.anodes];
         algn.ashifts = cellfunc(@(x) x - 1, obj.imagePositions);
         algn.aposition = obj.absoluteImagePositions;
         algn.aposition = algn.aposition{1};
         
         algn.irawdatasize = algn.absoluteSize;         
         algn.initializeCellDataSizeAndFormatFromRaw;
      end

      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % rois

      function [obj, roi] = reduceToROI(obj, roi)   
         [ids, shids] = obj.nodesFromROI(roi.boundingbox);
               
         obj.nodes = ids;
         obj.reducePairs;

         sh = cell2mat(obj.ashifts(shids));
         sh = min(sh, [], 1);
         roi.shift(-sh);
         
         obj.optimizePariwiseShifts;
      end
      
      
      function [nds, ids] = nodesFromROI(obj, roi)
         %
         % [nds, ids] = nodesFromROI(obj, roi)
         %
         % description:
         %    identifies the nodes that contribute to the roi
         %
         % input:
         %    roi    region of interest
         %
         % output:
         %   ids     nodes
         %   shids   ids of the nodes
         
         if isempty(obj.dataSize)
            error('%s: images not aligned !', class(obj));
         end
         
         ids = roiToImageIndex(obj.imagePositions, obj.imageSizes, roi);
         nds = obj.nodes;
         nds = nds(ids);
      end

     
%       function r = roi(obj, varargin)
%          %
%          % r = roi(obj)
%          %
%          % description:
%          %    generate outline of the component 
%       end       

      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % previews

      function p = preview(obj, varargin)
         for i = 1:length(obj)
            p = obj(i).asource.preview(obj(i).anodes);
         end
      end

      function p = previewStiched(obj, varargin)

         if length(obj) == 1  
            p = stitchPreview(obj, varargin);  
         else
            nobj = length(obj);
            imgs = cell(1, nobj);
            for i = 1:nobj
               imgs{i} = obj(i).previewStiched(varargin);
            end
            
            scale = obj(1).apreviewscale;
            shifts = obj.absolutePosition;
            shifts = cellfunc(@(x) round(x .* scale), shifts);
            %var2char(shifts)
            
            p = stitchPreview(imgs, 'shifts', shifts, 'scale', 1);
         end
      end
      
      function obj = clearPreview(obj)
         obj.asource.clearPreview();
         clearPreview@ImageSource(obj);
      end
      
      function obj = setPreviewScale(obj, scale)
         obj.asource.setPreviewScale(scale);
         setPreviewScale@ImageSource(obj, scale);
      end
      
%       function s = previewScale(obj)
%          s = obj.asource.previewScale();
%       end


      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % visualization

      function plotAlignedImages(obj)
         %
         % plotAlignedImages(obj)
         %
         % description:
         %    visualizes the alignment using plotAlignedImages
         %
         % See also: plot
         
         plotAlignedImages(obj.images, obj.imageShifts);
      end
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % information 
      
      function istr = infoString(obj)
         istr = infoString@ImageSource(obj, 'Aligned');
         istr = [istr, '\nsource:     ', class(obj.asource)];
         istr = [istr, '\nnodes:      ', var2char(obj.nNodes)];
         istr = [istr, '\npairs:      ', var2char(obj.nPairs)];
         istr = [istr, '\nalignerror: ', var2char(obj.alignmentError)];
      end

   end
   

   
   methods(Access = protected)
       % Override copyElement method:
      function cpObj = copyElement(obj)
         % copy of all four properties
         cpObj = copyElement@matlab.mixin.Copyable(obj);
         % Make a deep copy of the pairs object
         cpObj.apairs = copy(obj.apairs);
      end
   end
   
end