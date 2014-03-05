function seplabel = imlabelseparate(label)
%
% seplabel = imlabelseparate(label)
%
% description:
%     relabels the image by separating non touching regions into separat
%     of same into multiple label.
%
% input:
%     label       labeled image
% 
% output:
%     seplabel    relabled image, no separate regoins with same label
%
% note:
%     bwlabeln(label > 0) fails if two labels touch !
%
% See also: bwlabeln, bwlabel

labs = imlabel(label);
nextlabel = max(labs) + 1;
seplabel = label;

for l = labs
   cc = bwconncomp(label == l);
   for nl = 2:cc.NumObjects
      seplabel(cc.PixelIdxList{nl}) = nextlabel;
      nextlabel = nextlabel + 1;
   end
end

end
      
      
      