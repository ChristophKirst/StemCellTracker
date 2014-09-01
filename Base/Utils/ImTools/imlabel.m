function label = imlabel(labeledimage)
%
% label = imlabel(labeledimage)
%
% description:
%    return the list of lables in labeled image, excluding background = 0
%

label = unique(labeledimage(:))';
if ~isempty(label) && label(1) == 0
   label = label(2:end);
end
