function shifts = alignImages(imgs, varargin)
%
% shifts = alignImages(imgs, param)
%
% desciption:
%    infers shifts of images in cell array imgs to make them align
%    the fucntion is a wrapper that distributes the tasks to the indibvidual align and optimization methods
% 
% input:
%    imgs    images or image reading commands as cell array, if pairs parameter is not given the images are taken to 
%            be prealigned as in the cell outline, the cell coordinates are imgs{x,y,z}
%    param   (optional) parameter struct with entries:
%            .pairs           (optional) array struct with entries .from .to indicating the pairs (optional entry .orientation)
%            .method          'full'   = perform full alignment using key point detection and optimization using Hugin
%                             'global' = perform alignment given the individual pairwise alignments ('global')
%            .alignment       'Optimization', 'RMS', 'Correlation', 'Hugin' ('RMS')
%           
%            parameters for methods align2ImagesOnGridByXXX
%            other parameter as used by the various subroutines: alignImagesOnGird
%
% output: 
%    shifts  absolute shifts of images assuming [lower,left(,bottom)] of image in cell array starts a [0,0(,0)] 
%            in pixel units and pixel coordinates
%
% See also: alignImages, plotAlignedImages


param = parseParameter(varargin{:});

if ~iscell(imgs)
   error('alignImages: expecting cell array as input imgs');
end

adim = ndims1(imgs);
if adim > 3 || adim < 1
   error('alignImagesOnGrid: expect 1d - 3d grid of images, found %gd', adim);
end


methd = getParameter(param, 'method', 'global');

if strcmp(methd, 'full')
   shifts = alignImagesByHugin(imgs, param);
   
else
   
   pairs = getParameter(param, 'pairs', {});
 
   if isempty(pairs)
      shifts = alignImagesOnGrid(imgs, param);
   else

      if ~strcmp(methd, 'global')
         warning('alignImages: method %s not applicable with specified pairs: using global instead!', methd);
      end
      
      shifts = alignByGlobal(imgs, pairs, param);
   end
 
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% global shift optimization using predefined pairs

function shifts = alignByGlobal(imgs, pairs, param)

      alignnames = {'Optimization', 'Correlation', 'RMS', 'Hugin'};
      algn = getParameter(param, 'alignment', 'RMS');

      if ~any(strcmp(alignnames, algn))
         error('alignImages: alignment %s not in %s', align, var2char(alignnames));
      end

      fun = str2func(['align2ImagesBy' algn]);

      np = length(pairs);
      dim = ndims(imgs{1});
      sh(np) = struct('from', -1, 'to', -1, 'shift', zeros(1,dim));
      
      for p = 1:np
         pairs(p).shift =  fun(imgs{pair(p).from}, imgs{pair(2)}, param);
      end
 
      % globally optimize pairwise shifts

      [~, ic] = optimizePairwiseShifts(sh);
      ic = round(ic);
   
      % transform consistent shifts to shifts 
   
      shifts = num2cell(ic',2);
      shifts = zero2first(shifts);
      shifts = reshape(shifts, size(imgs));
end


% zero first to first shift
function shifts = zero2first(shifts)

   sh = shifts{1};
   for i = 1:numel(shifts)
      shifts{i} = shifts{i} - sh;
   end

end
