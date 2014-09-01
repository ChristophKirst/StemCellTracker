function ws = imlabelwatershed(image, label)
%
% ws = imlabelwatershed(image, label)
%
% description:
%    seed watershed with possible labeling
%
% input:
%   image    image to warershed
%   label    labele image
%
% output:
%   ws        watershedded image labeled according to the label input
%
% See also: watershed


labels = imlabel(label);
bd = impixelsurface(label);

% put seeds
lbl = label > 0;
for l = labels
   objb = (bd == l);
   if sum(bd(:)) > 1 % avoid seeds of single pixel to get lost
      lbl = lbl - objb;
   end   
   ll = ll + 1;
end

% water shed
image = mat2gray(image + 1 - lbl);
wsi = watershed(image);

%relabel
ws = zeros(size(image));

for l = labels(end:-1:1)
   obj = (label == l);
   idx = find(obj .* cast(wsi, class(label)), 1, 'first');
   wl = wsi(idx);   
   ws(wsi == wl) = l;
end

% correct for lost lables on boundaries

ws = ws + (wsi == 0) .* label;

end
