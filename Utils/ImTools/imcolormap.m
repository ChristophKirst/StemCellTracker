function varargout = imcolormap(varargin)
%
% colmap = imcolormap(varargin)
% colmap = imcolormap()
% colmap = imcolormap(ncols)
% colmap = imcolormap(colmapname)
% colmap = imcolormap(colmapname, ncols)
% colmap = imcolormap(..., perm)
%
% description: 
%     generates colormaps from a name, a number of labels or both
%
% input:
%     ncols       (optional) number of of colors (64)
%     colmapname  (optional) name of colormap ('default' = 'jet')
%     perm        (optional) permutation, either array of indices, 'random', 'shuffle' or 'none' ('none')
%
% output:
%     colmap      (optional) the color map as array, if nargout == 0, colormap(colmap) is called
%
% note: 
%     'shuffle' generates a random but fixed permuation, 'random' generates a varying permuation
%
% See also: colormap, imcolorize, imgray2color

colmapname = 'default';
ncol = 64;
perm = 'none';
colmap = [];

permnames = {'none', 'random', 'shuffle'};
usercolmapnames = {'igray', 'inversegray', 'default', 'y', 'yellow', 'm', 'magenta', 'c', 'cyan', 'r', 'red', 'g', 'green', 'b', 'blue', 'w', 'white', 'k', 'black'};
if nargin >=1
   var1 = varargin{1};  
   if isnumeric(var1)
      if size(var1,2) == 1
         ncol = var1;
      else
         colmap = var1;
         if ~ischar(colmap)
            ncol = size(colmap,1);
         end
         if size(colmap,2) ~= 3
            error('imcolormap: colormap needs 3 columns');
         end
      end
   elseif ischar(var1)
      if any(strcmp(permnames, var1))
         perm = var1;
      else
         colmapname = var1;
      end
   end
end


if nargin >=2
   var2 = varargin{2};
   if isnumeric(var2)
      ncol = var2;
   elseif ischar(var2)
      if any(strcmp(permnames, var2))
         perm = var2;
      else
         colmapname = var2;
      end
   end
end

if nargin == 3
   perm = varargin{3};
end


if isempty(colmap)
   if any(strcmp(usercolmapnames, colmapname))
      switch colmapname
         case {'igray', 'inversegray'}
            colmap = gray(ncol);
            colmap = flip(colmap, 1);
            % space for other nice color maps / colormaps that map to ncol etc... 
         case {'y', 'yellow', 'm', 'magenta', 'c', 'cyan', 'r', 'red', 'g', 'green', 'b', 'blue', 'w', 'white', 'k', 'black'}
            rgb = double(imcolorspec2rgb(colmapname));
            colmap = (double(0:(ncol-1)) / double(ncol-1))' * rgb;
         otherwise
            colmap = jet(ncol);   
      end
   else
      colmapname = str2func(colmapname);
      colmap = colmapname(ncol);
   end
end

%reduce or enhance range
ncm = size(colmap, 1);
if size(colmap, 1) ~= ncol
   ids = 0:ncol-1;
   colmap = colmap(min(round(ids/(ncol-1) * (ncm-1)+1),ncm), :);
end
     
% shuffle the data if requested
switch perm
   case 'random'     
      colmap = colmap(randperm(ncol), :);
   case 'shuffle'
      savedState  = RandStream.setGlobalStream(RandStream('mt19937ar', 'Seed', 1));
      colmap = colmap(randperm(ncol), :);
      RandStream.setGlobalStream(savedState); 
end
      

if nargout < 1
   colormap(colmap)
else
   varargout{1} = colmap;
end
   
end
   
