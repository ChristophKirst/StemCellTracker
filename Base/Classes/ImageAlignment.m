classdef ImageAlignment < matlab.mixin.Copyable
   %
   % ImageAlignment class representing alignment/tiling information for use with ImageSourceTiling class
   %

   properties
      ipairs = [];     % pairs of overlapping images (array of AlignmentPair classes)
      inodes = [];     % ids of images used for this alignment

      iasize = [];     % absolute image size of aligned images
   end
  
   methods

      function obj = ImageAlignment(varargin)
         %
         % ImageAlignment()
         % ImageAlignment(cellarray)
         % ImageAlignment(struct)
         % ImageAlignment(..., fieldname, fieldvalue, ...)
         %
         if nargin == 1 
            if isa(varargin{1}, 'ImageAlignment') %% copy constructor
               obj = copy(varargin{1});
            elseif isa(varargin{1}, 'AlignmentPair')
               obj.fromAlignmentPair(varargin{1});
            elseif isa(varargin{1}, 'cell')
               obj.fromCell(varargin{1});
            elseif isa(varargin{1}, 'struct')
               obj.fromStruct(varargin{1});
            elseif isnumeric(varargin{1})
               obj.fromTilingSize(varargin{1});
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
         obj.fromTilingSize(si);
      end
      
      
      function obj = fromStruct(obj, st)
         %
         % obj = fromStruct(obj, st)
         %
         % descritpion:
         %   initializes class from a struct st 
         %
         % input:
         %   st  struct with .from, .to and additional information such as shift, orientation, quality
         %
         
         obj.ipairs = AlignmentPair(st);
         obj.inodes = obj.uniquePairs();
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

         obj.ipairs = AlignmentPair(ap);  % deep copy
         obj.inodes = obj.uniquePairs();
      end
      
      
      function obj = fromTilingSize(obj, si) % assumes image ids arranged in array by reshape(1:pord(si), si)
         %
         % obj = fromTilingSize(obj, si)
         %
         
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
         
         obj.ipairs = p;
         obj.inodes = 1:prod(si);
      end
      
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% connectivity
      
      function obj = overlapQuality(obj, source, varargin)
         %
         % obj = overlapQuality(obj)
         %
         % description:
         %    calculates operlap quality of the images
         %
         % See also: overlapQuality, overlapStatisticsImagePair
         
         for p = 1:obj.npairs
            pp = obj.ipairs(p).copy();
            pp.from = source.getTile(pp.from);
            pp.to   = source.getTile(pp.to); 
            obj.ipairs(p).quality = overlapQuality(pp, varargin{:});
         end
      end
      
 
      function obj = removeLowQualityPairs(obj, thq)
         %
         % obj = removeLowQualityPairs(obj, thq)
         %
         % descritpion:
         %   removes pairs whose quality is less than thq
         %
         % input:
         %   thq   quality threshold 
         
         q = [obj.ipairs.quality];
         obj.ipairs = obj.ipairs(q >= thq); 
      end
      
      function obj = reducePairs(obj)
         %
         % obj = reduce(obj)
         %
         % description: 
         %     removes all paris not connecting two nodes in obj.nodes
         %

         iids = obj.inodes;
          
         id  = [obj.ipairs.from];
         ids = arrayfun(@(x) any(x==iids), id);
                  
         id  = [obj.ipairs.to];
         ids = and(ids, arrayfun(@(x) any(x==iids), id));
         
         obj.ipairs = obj.ipairs(ids);
      end
         

      function obj = reduceImages(obj)
         %
         % obj = reduceImages(obj)
         %
         % description: 
         %    relables nodes to 1:nnodes in nodes and pairs and orders/drops images accordingly
         %

         % make sure alignment is consistent
         if ~isempty(setdiff(obj.uniquePairs, obj.inodes))
            error('Alignment: there are image ids in pairs tthat do not appear in nodes: %s', var2char(setdiff(obj.uniquePairs, obj.inodes)));
         end
         
         iids = obj.inodes;
         nids = 1:length(iids);
         
         idconv = zeros(1,max(iids));
         idconv(iids) = nids;
         
         if ~isempty(obj.ipairs)
            pi = [obj.ipairs.from];
            pi = num2cell(idconv(pi));
            [obj.ipairs.from] = pi{:};
         
            pi = [obj.ipairs.to];
            pi = num2cell(idconv(pi));
            [obj.ipairs.to] = pi{:};
         end

         obj.inodes = nids;
         
         obj.iasize = [];
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


      
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% dims / ids
      
      function d = dim(obj)
         %
         % d = dim(obj)
         %
         % descritpion:
         %   number of single tile data dimensions
         %

         if ~isempty(obj.pairs)
            d = obj.ipairs(1).dim();
         else
            d = 2;
         end
      end
      
      function n = npairs(obj)
         %
         % n = npairs(obj)
         %
         n = length(obj.ipairs);
      end
      
      function n = nnodes(obj)
         %
         % n = nnodes(obj)
         %
         n = length(obj.inodes);
      end

      function id = uniquePairs(obj)
         %
         % id = uniquePairs(obj)
         %
         % descritpion:
         %   returns all ids used in pairs
         
         id = unique([[obj.ipairs.from], [obj.ipairs.to]]);
      end
            
      
      
      
      
      
      
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% alignment routines
      
      function shifts = imageShifts(obj, varargin)  % image shifts from pairwise shifts
         %
         % shifts = imageShifts()
         %
         % description
         %      image shifts from pairwise shifts

         if obj.nnodes == 1
            shifts = {zeros(1,obj.dim)};
            return
         end
         
         if nargin < 2
            anchor = 1;
         else
            anchor = varargin{1};
         end
         
         % make shifts consistent by optimizing them
         [~, ic] = optimizePairwiseShifts(obj.ipairs);
         ic = round(ic);
   
         % transform consistent shifts to shifts 
         shifts = num2cell(ic',2);
         shifts = shifts(obj.inodes);
         
         sh = shifts{anchor};
         for i = 1:numel(shifts)
            shifts{i} = shifts{i} - sh;
         end
      end
      
      
      function obj = alignPairsFromShifts(obj, ishifts)
         %
         % obj = alignPairsFromShifts(obj, ishifts)
         %
         % description:
         %    sets the pairwise shifts form image shifts
         %
         
         obj.ipairs = alignPairsFromShifts(obj.ipairs, ishifts);   
         
      end
      
      function obj = absoluteShiftsAndSize(obj, source)
         %
         % obj = absoluteImageShiftsAndSize(obj, source)
         %
         % description:
         %    calculates absolute size and shifts
         %
         % See also: absoluteShiftsAndSize
         
         [ashifts, as] = absoluteShiftsAndSize(obj.imageShifts, source.getTileSizes);
         
         obj.iasize = as;
         obj.alignPairsFromShifts(ashifts);
      end
      
      
      function obj = optimizePairwiseShifts(obj)
         %
         % obj = optimizePairwiseShifts(obj)
         %
         % description:
         %    globally optimizes pairwise shifts
         
         obj.ipairs = optimizePairwiseShifts(obj.pairs);
      end
      
      function obj = makeShiftsConsistent(obj)
         %
         % obj = makeShiftsConsistent(obj)
         %
         % description:
         %    makes shifts mutually consistent (i.e. paths in the grid commute)
         
         obj.ipairs = optimizePairwiseShifts(obj.pairs);
      end

      function obj = alignPairs(obj, source, varargin)
         %
         % alignPairs(obj, varargin)
         %
         % descritpion:
         %   alignes the individual paris of images
         %
         % input:
         %   param  parameter as for alignImagePair
         
         if isempty(source)
            error('alignment: alignParis: image suorce empty, canot align data!');
         end
         
         for i = 1:obj.npairs
            p = AlignmentPair(obj.ipairs(i));
            p.from = source.getTile(p.from);
            p.to   = source.getTile(p.to);
            
            p = alignImagePair(p, varargin{:});
           
            obj.ipairs(i).shift   = p.shift;
            obj.ipairs(i).quality = p.quality;
         end
      end
      

      function obj = align(obj, source, varargin)
         %
         % obj = align(obj, varargin)
         %
         % description:
         %    aligns images and sets new shifts
         %
         % See also: alignImages
         
         param = parseParameter(varargin);
         serial = getParameter(param, 'serial', false);
         
         if serial
            imgs = source;
         else
            imgs = source.getTiles;
         end

         [~, pairs] = alignImages(imgs, 'pairs', obj.ipairs, varargin{:});
         
         obj.ipairs = pairs;
         
         obj.absoluteShiftsAndSize(source);
      end
   
      function st = stitch(obj, source, varargin)
         %
         % st = stitch(obj, source, param)
         %
         % description
         %     stitches images using the alignment information and source
         
         st  = stitchImages(source.getTiles, obj.imageShifts, varargin{:});  
      end
      

      
      function plotAlignedImages(obj, source)
         %
         % plotAlignedImages(obj)
         %
         % description:
         %    visualizes the alignment using plotAlignedImages
         %
         % See also: plotAlignedImages
                  
         plotAlignedImages(source.getTiles, obj.imageShifts);
      end
      
   end
   
   
   methods(Access = protected)
       % Override copyElement method:
      function cpObj = copyElement(obj)
         % Make a shallow copy of all four properties
         cpObj = copyElement@matlab.mixin.Copyable(obj);
         % Make a deep copy of the pairs object
         cpObj.pairs = copy(obj.pairs);
      end
   end
   
   methods(Access = private)
   
   
   end
   
end