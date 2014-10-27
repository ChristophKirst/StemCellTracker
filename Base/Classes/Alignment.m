classdef Alignment < ImageSource
   %
   % Alignment class representing alignment/tiling data
   %

   properties
      apairs = [];       % pairs of overlapping images (array of AlignmentPair classes)
      anodes = [];       % ids of images used for this alignment
      
      asource = [];      % (optional) image source. assumed to have routines cell(id), cellSize(id) and data(id)

      aorigin = [0,0];   % absolute position of the first node
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
         if nargin == 1 
            if isa(varargin{1}, 'Alignment') %% copy constructor
               obj = copy(varargin{1});
            elseif isa(varargin{1}, 'AlignmentPair')
               obj.fromAlignmentPair(varargin{1});
            elseif isa(varargin{1}, 'cell')
               obj.fromCell(varargin{1});
            elseif isa(varargin{1}, 'struct')
               obj.fromStruct(varargin{1});
            elseif isa(varargin{1}, 'ImageSource')
               obj.fromImageSource(varargin{1});
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
         obj = obj.fromCellSize(size(ca));
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
      
      function obj = fromImageSource(obj, is)
         %
         % obj = fromImageSource(obj, is)
         %
         % descritpion:
         %   intitializes class from an ImageSource class
         %
         % input:
         %   is  ImageSource class
         
         obj.asource = is;
         obj.fromCellSize(is.cellSize);
      end
      
      
      function obj = initializeFromPairs(obj)
         % intializes the data size and format info using the info in pairs
         
         % TODO
         [~, obj.irawdatasize] = obj.absoluteShiftsAndSize;
         frmt = 'XYZ';
         frmt = frmt(1:length(obj.irawdatasize));
         obj.irawdataformat = frmt;
         
         obj.irawcellsize = 1;
         obj.irawcellformat = 'S';
         
         obj.initializeCellDataSizeAndFormatFromRaw;
      end
         

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % basics 
      
      function n = nodes(obj)
         n = obj.anodes;
      end
      
      function p = pairs(obj)
         p = obj.apairs;
      end
      
      function s = source(obj)
         s = obj.asource;
      end
      
      function o = origin(obj)
         o = obj.aorigin;
      end


      function d = dim(obj)
         %
         % d = dim(obj)
         %
         % descritpion:
         %   number of dimensions of the alignment shifts
         %

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
         n = length(obj.apairs);
      end
      
      function n = nNodes(obj)
         %
         % n = nnodes(obj)
         %
         n = length(obj.anodes);
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
         is = obj.asource.dataSize(varargin{:});
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
         obj = obj.cleatCache;
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

            % we dont want to stitch a huge image for a small roi -> find relevant ids and pairs first
            
            [ids, shids] = obj.nodeIndexFromROI(roi);

            ash = obj.imageShifts;
            ash = ash(shids);
            tsi =obj.imageSizes;
            tsi= tsi(shids);

            st  = stitchImages(obj.sourceCell(ids), ash, varargin{:});
            
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
         %   give an array of alignments remove the ones in the node list
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
         
         obj.apairs = optimizePairwiseShifts(obj.pairs);
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
         
         obj.initializeFromPairs;
      end

      function st = stitch(obj, varargin)
         %
         % st = stitch(obj, source nsubalgnnsubalgn , param)
         %
         % description
         %     stitches images using the alignment information and source
         
         st  = stitchImages(obj.images, obj.imagePositions, varargin{:});
         
         if obj.irawcache
            obj.fromData(st);
         end
      end
      
      function st = stitchArray(obj, varargin)
         % stitches an array of alignement classes using thier different origins
         st  = stitchImages([obj.images], [obj.imagePositions], varargin{:});
      end
      
      
      function obj = setOrigin(obj, o)
         obj.aorigin = o;
      end
      
      function obj = setShifts(obj, s)
         obj.ashifts = s;
      end
 
      %%% routines querrying exsiting results
      
      function oq = overlapQuality(obj)
         oq = [obj.apairs.quality];
      end
      
      
      function [ashifts, as] = absoluteShiftsAndSize(obj)
         %
         % obj = absoluteShiftsAndSize(obj)
         %
         % description:
         %    calculates absolute size and shifts
         %
         % See also: absoluteShiftsAndSize
         
         [ashifts, as] = absoluteShiftsAndSize(obj.imageShifts, obj.imageSizes);
         %obj.pairs = alignPairsFromShifts(obj.pairs, ashifts, obj.nodes);
      end

      function isizes = imageSizes(obj)
         isizes = repmat({obj.sourceDataSize}, 1, obj.nNodes); 
      end

      function imgs = images(obj)
         imgs = obj.asource.cell(obj.nodes);
      end

      function ipos = imagePositions(obj)
         %
         % ipos = imagePositions()
         %
         % description:
         %    returns positions of the images given the origin and shifts
         %
         % See also: optimizePairwiseShifts

         if obj.nNodes == 1
            ipos = {ones(1,obj.dim)};
            return
         end
         
         [~, ic] = optimizePairwiseShifts(obj.apairs);
         ic = round(ic);
   
         % transform consistent shifts to shifts 
         ipos = num2cell(ic',2);
         %var2char(shifts);
         
         ipos = cellfunc(@(x) x + obj.aorigin, ipos);
      end

      function shifts = imageShifts(obj)
         %
         % shifts = imageShifts()
         %
         % description:
         %    returns absolute shifts of each image determined from pairwise shifts
         %
         % See also: optimizePairwiseShifts

         if obj.nNodes == 1
            shifts = {zeros(1,obj.dim)};
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
         
         comp = connectedAlignments(obj, varargin{:});   
      end
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % utils
      
      function [obj, roi] = reduceToROI(obj, roi)   
         [ids, shids] = obj.roi2tileids(roi.boundingbox);
      
         sh = cell2mat(obj.ashifts(shids));
         sh = min(sh, [], 1);
         
         obj.nodes = ids;
         obj.reducePairs;

         roi.shift(-sh);
         
         obj.absoluteShiftsAndSize;
      end
      
            
      function [ids, shids] = nodeIndexFromROI(obj, roi)
         %
         % [ids, shids] = roi2tileids(obj, roi)
         %
         % description:
         %    identifies the ids of the tiles that contribute to the roi
         %
         % input:
         %    roi    regoin of interest
         %
         % output:
         %   ids     node ids
         %   shids   id of the nodes
         
         if isempty(obj.dataSize)
            error('%s: images not aligned !', class(obj));
         end
         
         shids = roiToImageIndex(obj.imagePositions, obj.imageSizes, roi);
         nds = obj.nodes;
         ids = nds(shids);
      end

      function aerror = alignmentError(obj)
          if obj.nPairs == 0
              aerror = 0;
          else
            aerror = sum([obj.apairs.aerror]) / obj.nPairs;
          end
      end

      
      function r = roi(obj, varargin)
         %
         % r = roi(obj)
         %
         % description:
         %    generate outline of the component
         %    if images are not aligned assume grid structure and max
         %    alignment
         
         
         
         
      end
      

      function ov = overlaps(obj, varargin)
         
         
         
      end
      
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % visualization
   
      
      function d = plotAlignedPreview(obj, varargin)
         d2 = alignCellPreview(obj, varargin);

         if nargout < 1
            implot(d2)
         else
            d = d2;
            implot(d)
         end
      end
  
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
      
      
      function plot(obj, varargin)
         %
         % plot(obj)
         %
         % description:
         %    plot the stitched image
         %
         % See also: plotAlignedImages
                  
         implot(obj.stitch(varargin{:}));
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