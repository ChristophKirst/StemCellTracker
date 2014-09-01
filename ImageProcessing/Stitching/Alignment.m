classdef Alignment < matlab.mixin.Copyable
   %
   % Alignment class representing alignment/tiling data
   %

   properties
      pairs = [];       % pairs of overlapping images (array of AlignmentPair classes)
      sizes  = {};      % image sizes as cell array
      asize = [];       % absolute image size of aligned images
      
      images = {};      % (optional) image data
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

      function d = dim(obj)
         %
         % d = dim(obj)
         %
         % descritpion:
         %   number of alignment dimensions
         %
         
         if ~isemtpy(obj.pairs)
            d = obj.pairs(1).dim();
         else
            d = 2;
         end
      end
      
      function n = npairs(obj)
         n = length(obj.pairs);
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
         obj.images = ca(:);
         obj.sizes  = cellfunc(@size, obj.images);
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
      end
      
      
      function sub = subAlignment(obj, pids)
         %
         % sub = subAlignment(obj, pids)
         %
         % descritpion:
         %   creates sub Alignment from pairs with ids pids
         %
         % input:
         %   pids  ids of pairs to use in sub
              
         sub = Alignment(obj.pairs(pids));
         
         sub.sizes = obj.sizes;
         sub.images = obj.images;
      end
      
      
      function obj = relabel(obj)
         %
         % obj = relabel(obj)
         %
         % description: 
         %     relables the nodes to the nodes used in pairs and reduces the images and sizes arrays accordingly
         %     resets asize
         
         iids = unique([[obj.pairs.from], [obj.pairs.to]]);
         nids = 1:length(iids);
         
         idconv = zeros(1,max(iids));
         idconv(iids) = nids;
         
         pi = [obj.pairs.from];
         pi = num2cell(idconv(pi));
         [obj.pairs.from] = pi{:};
         
         pi = [obj.pairs.to];
         pi = num2cell(idconv(pi));
         [obj.pairs.to] = pi{:};
         
         if ~isempty(obj.sizes)
            obj.sizes = obj.sizes(iids);  
         end
                  
         if ~isempty(obj.images)
            obj.images = obj.images(iids);  
         end
         
         obj.asize = [];
      end
      
      
      function obj = imageSizes(obj)
         %
         % obj = imageSizes(obj)
         %
         % descritpion:
         %   calcualtes the image sizes

         obj.sizes = cellfunc(@size, obj.images);
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
         
         if isempty(obj.images)
            error('alignment: alignParis: no image data to align images, set property images');
         end
         
         for i = 1:obj.npairs
            p = AlignmentPair(obj.pairs(i));
            p.from = obj.images{p.from};
            p.to   = obj.images{p.to};
            
            p = alignImagePair(p, varargin{:});
           
            obj.pairs(i).shift   = p.shift;
            obj.pairs(i).quality = p.quality;
         end
      end


      function comp = connectedAlignments(obj, varargin)
         %
         % comp = connectedAlignments(obj, varargin)
         %
         % descritpion:
         %   given quality factors of piars retunrs cell array of Alignemnts with connected images
         %
         % input:
         %   param   paramter as in connectedAlignments
         %
         % See also: conenctedAlignments
         
         comp = connectedAlignments(obj);
         
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
            stats = overlapStatisticsImagePair(obj.pairs(p), varargin{:});
            obj.pairs(p).quality = overlapQuality(stats, varargin{:});
         end
      end
      
      function obj = absoluteShiftsAndSize(obj)
         %
         % obj = absoluteShiftsAndSize(obj)
         %
         % description:
         %    calculates absolute size and shifts
         %
         % See also: absoluteShiftsAndSize
         
         [ashifts, as] = absoluteShiftsAndSize({obj.pairs.shifts}, obj.sizes);
         
         obj.asize = as;
         [obj.pairs.shifts] = ashifts{:};
      end
      
     
      function obj = optimizePairwiseShifts(obj)
         %
         % obj = optimizePairwiseShifts(obj)
         %
         % description:
         %    globally optimizes pairwise shifts by LMS
         %
         % See also: optimizePairwiseShifts
         
         obj.pairs = optimizePairwiseShifts(obj.pairs);
         
      end
      
      
      function obj = align(obj, varargin)
         %
         % obj = align(obj, varargin)
         %
         % description:
         %    aligns images and sets new shifts
         %
         % See also: alignImages
         
         shifts = alignImages(obj.images, 'pairs', obj.pairs, varargin{:});
         
         [obj.pairs.shift] = shifts{:};
      end
      

      function plot(obj)
         %
         % plot(obj)
         %
         % description:
         %    visualizes the alignment using plotAlignedImages
         %
         % See also: plotAlignedImages
         
         plotAlignedImages(obj.images, reshape({obj.pairs.shifts}, size(obj.images)));
      end
      
   end
      
end