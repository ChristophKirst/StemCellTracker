function [shifts, pairs] = alignImages(imgs, varargin)
%
% [shifts, pairs] = alignImages(imgs, param)
%
% desciption:
%    infers shifts of images in cell array imgs to make them align
%    the fucntion is a wrapper that distributes the tasks to the indibvidual align and optimization methods
% 
% input:
%    imgs    images or ImageSource class, if pairs parameter is not given the images are taken to 
%            be prealigned as in the cell outline, the cell coordinates are imgs{x,y,z}
%    param   (optional) parameter struct with entries:
%            .pairs           (optional) array struct with entries .from .to indicating the pairs (optional entry .orientation)
%            .method          'full'   = perform full alignment using key point detection and optimization using Hugin
%                             'global' = perform alignment given the individual pairwise alignments ('global')
%            .alignment       'Optimization', 'RMS', 'Correlation', 'Hugin' ('RMS')
%           
%            parameters used by align2ImagesOnGridByXXX
%
% output: 
%    shifts  absolute shifts of images assuming [lower,left(,bottom)] of image in cell array starts a [0,0(,0)] 
%            in pixel units and pixel coordinates
%    pairs   the aligned pairs
%
% See also: Alignment, align2ImagesOnGrid, plotAlignedImages


param = parseParameter(varargin{:});

istiling = false;
if iscell(imgs)
   adim = ndims1(imgs);
elseif isa(imgs, 'ImageSource')
   istiling = true;
   adim = length(imgs.cellsize);   
else
   error('alignImages: expecting cell array or ImageSource as input imgs');
end

if adim > 3 || adim < 1
   error('alignImagesOnGrid: expect 1d - 3d grid of images, found %gd', adim);
end


methd = getParameter(param, 'method', 'global');

if strcmp(methd, 'full')
   
   if istiling
      imgs = imgs.getTiles();
   end
   shifts = alignImagesByHugin(imgs, param);
   
   if nargout > 1
      pairs = alignPairsFromShiftsOnGird(shifts); % TODO: strictly speaking this needs to be more general not on gird !
   end
   
else
   
   pairs = getParameter(param, 'pairs', {});
 
   if isempty(pairs)
      if istiling
         imgs = imgs.getTiles();
      end
      shifts = alignImagesOnGrid(imgs, param);
      
      if nargout > 1
         pairs = alignPairsFromShiftsOnGird(shifts);
      end
      
   else

      if ~strcmp(methd, 'global')
         warning('alignImages: method %s not applicable with specified pairs: using global instead!', methd);
      end
      
      [shifts, pairs] = alignByGlobal(imgs, pairs, istiling, param);
   end
 
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% global shift optimization using predefined pairs

function [shifts, pairs] = alignByGlobal(imgs, pairs, istiling, param)

      alignnames = {'Optimization', 'Correlation', 'RMS', 'Hugin'};
      algn = getParameter(param, 'alignment', 'RMS');

      if ~any(strcmp(alignnames, algn))
         error('alignImages: alignment %s not in %s', align, var2char(alignnames));
      end



      np = length(pairs);
  
      if isentry(pairs, 'orientation') && ~isempty(pairs(1).orientation)  %% assuming either all pairs have oriantion or none
         isorient = true;
         fun = str2func(['align2ImagesOnGridBy' algn]);
      else
         isorient = false;
         fun = str2func(['align2ImagesBy' algn]);
      end  
      %isorient

      if istiling
         
         if isorient
            for p = 1:np
               switch pairs(p).orientation
                  case 1
                     aimgs = {imgs.getTile(pairs(p).from); imgs.getTile(pairs(p).to)};
                  case 2
                     aimgs = {imgs.getTile(pairs(p).from), imgs.getTile(pairs(p).to)};
                  case 3
                     c{1,1,1} = imgs.getTile(pairs(p).from); c{1,1,2} = imgs.getTile(pairs(p).to);
                     aimgs = c;
                  otherwise
                     error('alignImages: image pair %g does not have valid orientaiton!', p);
               end

               pairs(p).shift = fun(aimgs, param);
            end
            
         else
         
            for p = 1:np
               pairs(p).shift =  fun(imgs.getTile(pairs(p).from), imgs.getTile(pairs(p).to), param);
            end
         end
         
      else  % image cell array input 
         
         if isorient
            for p = 1:np
               switch pairs(p).orientation
                  case 1
                     aimgs = {imgs{pairs(p).from}; imgs{pairs(p).to}};
                  case 2
                     aimgs = {imgs{pairs(p).from}, imgs{pairs(p).to}};
                  case 3
                     c{1,1,1} = imgs{pairs(p).from}; c{1,1,2} = imgs{pairs(p).to};
                     aimgs = c;
                  otherwise
                     error('alignImages: image pair %g does not have valid orientaiton!', p);
               end

               pairs(p).shift =  fun(aimgs, param);
            end
         else

            for p = 1:np
               pairs(p).shift =  fun(imgs{pairs(p).from}, imgs{pairs(p).to}, param);
            end
         end
      end
 
      % globally optimize pairwise shifts

      [~, ic] = optimizePairwiseShifts(pairs);
      ic = round(ic);
   
      % transform consistent shifts to shifts 
   
      shifts = num2cell(ic',2);
      shifts = zero2first(shifts);
      
      
      if ~istiling
         si = size(imgs);
         %size(shifts)
         shifts = reshape(shifts, si);
      end

      %var2char(shifts)
      pairs = alignPairsFromShifts(pairs, shifts);
      
end


% zero first to first shift
function shifts = zero2first(shifts)

   sh = shifts{1};
   for i = 1:numel(shifts)
      shifts{i} = shifts{i} - sh;
   end

end
