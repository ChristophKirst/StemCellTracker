function implotcolormap(cmap, siz)
%
% implotcolormap(cmap, siz)
%
% description:
%     plots the color map in a array
%
% input:
%    cmpa   color map
%    siz    (optional) size of array to plot ([] = square, 'horizontal' = 'row', 'vertical = 'colum')
%
% See also: colormap

if nargin < 2
   siz = [];
end
if strcmp(siz, 'automatic')
   siz = [];
end
if strcmp(siz, 'row')
   siz = 'horizontal';
end
if strcmp(siz, 'column')
   siz = 'vertical';
end

n = size(cmap,1);
dat = 1:n;

if isempty(siz)
   w = floor(sqrt(n));
   h = ceil(n/w);
   while h*w < n
      h = h + 1;
   end
   nn = h*w;
   
   dat(n+1:nn) = 1;
   siz = [h, w];
elseif ischar(siz)
   if strcmp(siz, 'vertical')
      siz = [1,n];
   elseif strcmp(siz, 'horizontal')
      siz = [n,1];
   else
      error('plotcolormap: size spec string invalid.'); 
   end
else
   error('plotcolormap: size spec invalid.'); 
end
      
cdat = cmap(dat, :);
cdat = cat(3, reshape(cdat(:,1), siz)', reshape(cdat(:,2), siz)', reshape(cdat(:,3), siz)');
imagesc('CData', cdat)
set(gca,'YDir','normal')

end






