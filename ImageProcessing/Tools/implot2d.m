function img = implot2d(img, varargin)
%
% implot2d(img, varargin)
%
% description:
%    implot plots image using pixel coordinates [p,q]
%
% input:
%    img      image to show
%    varargin color.scale as [cmin, cmax] and further options pasased to imshow
%
% See also: implot3d, imshow


%in case we extrated a z slice -> squeeze
img = squeeze(img);

%colormap
if nargin >= 2 
   if ismatrix(varargin{1}) && size(varargin{1},2) == 3
      cmap = varargin{1};
      vararg = varargin(2:end);
   else
      cmap = colormap;
      vararg = varargin;
   end
else
   cmap = colormap;
   vararg = varargin;
end

clim = [min(img(:)), max(img(:))];
for i = 1:2:length(vararg)
   if ischar(vararg{i}) && strcmpi(vararg{i}, 'color.scale')
      clim = vararg{i+1};
      vararg([i, i+1]) = []; 
   end
end
clim = double(clim);

%permute to match pq coordinates
switch ndims(img)
   case 2 % xy
      ncols = size(cmap, 1);  
      img = ceil( (double(img) - clim(1)) / (clim(2)-clim(1)) * ncols);
      img(img<1) = 1;
      img(img>ncols) = ncols;
      
      subimage(permute(img, [2 1 3]), cmap)
      
      %imshow(permute(img, [2 1 3]), cmap, 'DisplayRange', [0,max(img(:))], vararg{:});
      %freezecolormap(gca);
   case 3 % pqc
      imshow(permute(img, [2 1 3]), vararg{:});
end
      
axis on
xlabel('X'); ylabel('Y'); 
set(gca,'YDir','normal');


end


%matlab does not display image if one wants to rotate it
%view(-90,90)

% convert to rgb snippet
% 
%     ncols = size(cmap, 1);  
%     idx = ceil( (double(img) - clim(1)) / (clim(2)-clim(1)) * nColors);
%     idx(idx<1) = 1;
%     idx(idx>ncols) = ncols;
% 
%     imgrgb = zeros(siz);
%     for i = 1:3,
%         c = cmap(idx,i);
%         c = reshape(c,siz);
%         c(nanmask) = nancolor(i); %restore Nan (or nancolor if specified)
%         realcolor(:,:,i) = c;
%     end
