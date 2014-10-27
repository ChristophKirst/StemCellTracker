classdef AlignmentPair < matlab.mixin.Copyable
   %
   % AlignmentPair class representing two images with overlap to be aligned
   %
   
   properties
      from = [];         % first image id
      to = [];           % second image id
      orientation = 0;   % orientation of images to be aligned: 'lr'=1=left-right, 'du'=2=down-up 'bt'=3=bottom-top, ''=0=no primary direction
      shift = [0,0];     % shift between the images
      quality = 1;       % quality of overlap used to detect empty image overlaps
      aerror  = Inf;     % error detected when aligninng the image pair
   end
  
   methods

      function obj = AlignmentPair(varargin)
         %
         % AlignmentPair()
         % AlignmentPair(from, to, orientation)
         % AlignmentPair(from, to, shift)
         % AlignmentPair(cellarray)
         % AlignmentPair(struct)
         % AlignmentPair(..., fieldname, fieldvalue, ...)
         %
         if nargin == 1 
            if isa(varargin{1}, 'AlignmentPair') %% copy constructor
               obj = copy(varargin{1});
            elseif isa(varargin{1}, 'cell')
               obj.fromCell(varargin{1});
            elseif isa(varargin{1}, 'struct')
               obj.fromStruct(varargin{1});
            else
               error('%s: not valid arguments for constructor', class(obj));
            end
         elseif nargin == 2 && all(cellfun(@isnumeric, varargin(1:2)))
            obj.from = varargin{1};
            obj.to   = varargin{2};
         elseif nargin == 3 && all(cellfun(@isnumeric, varargin(1:2)))
            obj.from = varargin{1};
            obj.to   = varargin{2};
            if ischar(varargin{3})
               obj.orientation = varargin{3};  
            else
               obj.shift = varargin{3};  
            end
         else
            obj.fromParameter(varargin);
         end
         
         obj.initialize();
      end
      
      function initialize(obj)
         for i = 1:length(obj)
            obj(i).orientation = obj.orientation2number(obj(i).orientation);
         end
      end
      
      function obj = fromParameter(obj, varargin)
         obj = classFromParameter(obj, [], varargin);
         obj.initialize();
      end
      
      function obj = fromStruct(obj, st)
         %
         % obj = fromStruct(obj, st)
         %
         % descritpion:
         %   initializes AlignmentPair array form struct array st
         
         for i = 1:length(st)
            obj(i).from  = st(i).from;         
            obj(i).to    = st(i).to;        
         end
         
         fields = {'orientation', 'shift', 'quality'};
         for f = 1:length(fields)
            fn = fields{f};
            if isfield(st, fn)
               for i = 1:length(st)
                  obj(i).(fn) = st(i).(fn);
               end
            end
         end  
         
         obj.initialize();
      end

      function obj = fromCell(obj, ca)         
         %
         % obj = fromCell(obj, st)
         %
         % descritpion:
         %   initializes AlignmentPair form a pair of ids given as peraligned cell array
         
         si = size(ca);
         obj.from = ca{1};
         obj.to   = ca{2};
         obj.orientation = find(si == 2, 1);
      end
      
      function ca = toCell(obj)         
         %
         % ca = toCell(obj) 
         %
         % descritpion:
         %   converts the data into cell array representation
         
         if obj.orientation > 0
            pos = num2cell(ones(1,obj.dim));
            ca{pos{:}} = obj.from;
            pos{obj.orientation} = 2;
            ca{pos{:}} = obj.to;
         else
            ca = {obj.from, obj.to};
         end
      end
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % basics

      function d = dim(obj)
         %
         % d = dim(obj)
         %
         % descritpion:
         %   dimension of alignment problem = dimension of shift vector
         
         d = ndims(obj.shift);
      end
      
      function n = orientation2number(~, o)
         if ischar(o)
            switch o
               case 'lr'
                  n = 1;
               case 'du'
                  n = 2;
               case 'bt'
                  n = 3;
               case ''
                  n = 0;
               otherwise
                  error('AlignmentPair: orientation %s is not lr, du or bt or empty', o);
            end
         elseif isempty(o)
            n = 0;
         elseif isnumeric(o)
            n = o;
         else
            error('AlignmentPair: orientation is not lr, du or bt or empty: %s', var2char(o));
         end
      end

      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % alignment
      
      function obj = align(obj, varargin)
         %
         % obj = align(obj, varargin)
         %
         % description:
         %    aligns the two images using alignImagePair
         %
         % See also: alignImagePair
         
         obj = alignImagePair(obj, varargin{:});
      end
      
      
      function [ovl1, ovl2] = overlap(obj)
         %
         % [ovl1, ovl2] = overlap(obj)
         %
         % description:
         %    returns overlap region of the two images from and to using shift
         %
         % See also: overlap2AlignedImages
         
         [ovl1, ovl2] = overlap2AlignedImages(obj.from, obj.to, obj.shift);
      end
      
      function [ovl1, ovl2] = overlapOnGrid(obj, varargin)
         %
         % [ovl1, ovl2] = overlapOnGird(obj)
         %
         % description:
         %    returns overlap region of the two images using the overlap.max parameter and orientation
         %
         % See also: overlap2ImagesOnGird
         
         [ovl1, ovl2] = overlap2ImagesOnGrid(obj.toCell, varargin{:});       
      end
      
      
      
      
      function ovl = overlapSize(obj, imSizes)
         %
         % ovl = overlapSize(obj, imSizes)
         %
         % description:
         %    returns overlap sizes
         
         nobj = length(obj);
         ovl = zeros(nobj, length(obj(1).shift));
         
         for i = 1:nobj
            if nargin == 1
               imSizes = {size(obj(i).from), size(obj(i).to)};
            end

            sh = obj(i).shift;
            for j = 1:length(sh)
               if sh(j) >= 0
                  ovl(i, j) = max(0, imSizes{1}(j) - sh(j));
               else
                  ovl(i, j) = max(0, imSizes{2}(j) + sh(j));
               end
            end
         end
      end
      
      function ovl = overlapSizePrimary(obj, imSize)
         nobj = length(obj);
         ovl = zeros(1, nobj);
         for i = 1:nobj
            o = obj(i).orientation;
            if o > 0
               ovl(i) = imSize(o) - obj(i).shift(o);
            end
         end
         
      end


      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % visualization
      
      function plot(obj)
         %
         % plot(obj)
         %
         % description:
         %    visualizes the alignment using plot2AlignedImages
         %
         % See also: plot2AlignedImages
         
         plot2AlignedImages(obj.from, obj.to, obj.shift)
      end
      
   end
      
end