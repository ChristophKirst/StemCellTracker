classdef Alignment < matlab.mixin.Copyable
   %
   % Alignment class representing alignment/tiling data
   %

   properties
      pairs = [];      % pairs of overlapping images (array of AlignmentPair classes)
      nodes = [];      % ids of images used for this alignment
      
      isource = [];    % (optional) image data sources. assumed to have routines tile(id), tilesize(id)
   end
  
   methods

      function obj = Alignment(varargin)
         %
         % Alignment()
         % Alignment(cellarray)
         % Alignment(struct)
         % AlignmentPair(..., fieldname, fieldvalue, ...)
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
            elseif isa(varargin{1}, 'ImageSourceTiled')
               obj.fromImageSourceTiled(varargin{1});
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
         end
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
         
         obj.pairs = p;
         obj.nodes = 1:numel(ca);
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
         
         obj.pairs = AlignmentPair(st);
         obj.nodes = obj.pairIds();
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

         obj.pairs = AlignmentPair(ap);  % deep copy
         obj.nodes = obj.pairIds();
      end
      
      function obj = fromImageSourceTiled(obj, ist)
         %
         % obj = fromImageSourceTiled(obj, ist)
         %
         % descritpion:
         %   intitializes class from ImageSourceTiled class
         %
         % input:
         %   ist  ImagesourceTiled class
         
         obj.isource = ist;
         ts = ist.tileSizes();
         obj.fromCell(ts);
      end
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % basics 
      
      function d = dim(obj)
         %
         % d = dim(obj)
         %
         % descritpion:
         %   number of dimensions of the alignment shifts
         %

         if ~isempty(obj.pairs)
            d = obj.pairs(1).dim();
         else
            d = 2;
         end
      end
      
      function n = npairs(obj)
         %
         % n = npairs(obj)
         %
         n = length(obj.pairs);
      end
      
      function n = nnodes(obj)
         %
         % n = nnodes(obj)
         %
         n = length(obj.nodes);
      end

      function id = pairIds(obj)
         %
         % id = pairIds(obj)
         %
         % descritpion:
         %   returns all ids used in pairs
         
         id = unique([[obj.pairs.from], [obj.pairs.to]]);
      end      
      
      function ts = tilesizes(obj)
         ts = obj.isource.tilesizes;
         ts = ts(obj.nodes);
      end
      
      function tf = tileformat(obj)
         tf = obj.isource.tileformat;
      end
      
      
      function img = image(obj, id)
         nds = obj.nodes;
         img = obj.isource.tile(nds(id));
      end

      function is = imageSizes(obj, varargin)
         is = obj.isource.tilesizes;
         nds = obj.nodes;
         if nargin > 1
            is = is(nds(varargin{1}));
         else
            is = is(nds);
         end
      end
         
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % alignment routines
      
      
      % routines modifying the class
      
      function obj = removeLowQualityPairs(obj, thq)
         %
         % obj = removeLowQualityPairs(obj, thq)
         %
         % descritpion:
         %   removes pairs whose quality is less than thq
         %
         % input:
         %   thq   quality threshold 
         
         q = [obj.pairs.quality];
         obj.pairs = obj.pairs(q >= thq); 
      end
      
      function obj = reducePairs(obj)
         %
         % obj = reduce(obj)
         %
         % description: 
         %     removes all paris not connecting two nodes in obj.nodes
         %

         iids = obj.nodes;
          
         id  = [obj.pairs.from];
         ids = arrayfun(@(x) any(x==iids), id);
                  
         id  = [obj.pairs.to];
         ids = and(ids, arrayfun(@(x) any(x==iids), id));
         
         obj.pairs = obj.pairs(ids);
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

         
         nds = obj.nodes;
         n = obj.nnodes;
         rmn = zeros(1,n);
         for i = 1:n
            img = obj.isource.tile(nds(i));
            rmn(i) = max(img(:)) < thb;
         end
         
         if any(rmn)
            obj.nnodes = nds(~rmn);
            obj.reducePairs();
         end
      end
      
      function obj = overlapQuality(obj, varargin)
         %
         % obj = overlapQuality(obj) 
         %
         % description:
         %    calculates operlap quality of the images
         %
         % See also: overlapQuality, overlapStatisticsImagePair
         
         for p = 1:obj.npairs
            pp = obj.pairs(p).copy();
            pp.from = obj.isource.getTile(pp.from);
            pp.to   = obj.isource.getTile(pp.to); 
            obj.pairs(p).quality = overlapQuality(pp, varargin{:});
         end
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
         
         if isempty(obj.isource)
            error('alignment: alignParis: no image data to align images, set property isource');
         end
         
         for i = 1:obj.npairs
            p = AlignmentPair(obj.pairs(i));
            p.from = obj.isource.getTile(p.from);
            p.to   = obj.isource.getTile(p.to);
            
            p = alignImagePair(p, varargin{:});
           
            obj.pairs(i).shift   = p.shift;
            obj.pairs(i).quality = p.quality;
         end
      end
     
      function obj = optimizePairwiseShifts(obj)
         %
         % obj = optimizePairwiseShifts(obj)
         %
         % description:
         %    globally optimizes pairwise shifts 
         %
         % See also: optimizePairwiseShifts
         
         obj.pairs = optimizePairwiseShifts(obj.pairs);
      end


      % routines not modifying the class
      
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
 
      function [ashifts, as] = absoluteShiftsAndSize(obj)
         %
         % obj = absoluteShiftsAndSize(obj)
         %
         % description:
         %    calculates absolute size and shifts
         %
         % See also: absoluteShiftsAndSize
                
%          obj.nodes
%          obj.tilesizes
%          obj.imageShifts
         
         [ashifts, as] = absoluteShiftsAndSize(obj.imageShifts, obj.tilesizes);
         
         %obj.pairs = alignPairsFromShifts(obj.pairs, ashifts, obj.nodes);
      end
   
      function shifts = imageShifts(obj)
         %
         % shifts = imageShifts()
         %
         % description:
         %    returns absolute shifts of each image determined from pairwise shifts
         %
         % input:
         %    anchor   (optional) = 1
         %

         if obj.nnodes == 1
            shifts = {zeros(1,obj.dim)};
            return
         end
         
         [~, ic] = optimizePairwiseShifts(obj.pairs);
         ic = round(ic);
   
         % transform consistent shifts to shifts 
         shifts = num2cell(ic',2);
         %var2char(shifts);
         
         %todo optimize optimization to nodes that are really there and not max node number
         shifts = shifts(obj.nodes);
         %var2char(shifts);
         
         sh = shifts{1};
         for i = 1:numel(shifts)
            shifts{i} = shifts{i} - sh;
         end
      end
      
      
      % align and stitch

      function obj = align(obj, varargin)
         %
         % obj = align(obj, varargin)
         %
         % description:
         %    aligns images and sets new shifts
         %
         % See also: alignImages
         
         [~, prs] = alignImages(obj.isource, 'pairs', obj.pairs, varargin{:});
         obj.pairs = prs;
      end

      function st = stitch(obj, varargin)
         %
         % st = stitch(obj, source, param)
         %
         % description
         %     stitches images using the alignment information and source
         
         st  = stitchImages(obj.isource.getTiles(obj.nodes), obj.imageShifts, varargin{:});
      end
      

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % plotting
      
      function plotAlignedImages(obj)
         %
         % plotAlignedImages(obj)
         %
         % description:
         %    visualizes the alignment using plotAlignedImages
         %
         % See also: plot
         
         plotAlignedImages(obj.isource.tiles(obj.nodes), obj.imageShifts);
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
      
   end
   
   
   methods(Access = protected)
       % Override copyElement method:
      function cpObj = copyElement(obj)
         % copy of all four properties
         cpObj = copyElement@matlab.mixin.Copyable(obj);
         % Make a deep copy of the pairs object
         cpObj.pairs = copy(obj.pairs);
      end
   end
   
%    methods(Access = private)
%       
%    end
   
end