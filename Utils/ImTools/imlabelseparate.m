function seplabel = imlabelseparate(label)
%
% seplabel = imlabelseparate(label)
%
% description:
%     relabels the image by separating non touching regions with same label
%     into multiple label.
%
% input:
%     label       labeled image
% 
% output:
%     seplabel    relabled image, no separate regoins with same label
%
% note:
%     bwlabeln(label > 0) fails if two labels touch
%     bwconncomp scales badly with image size even for small objects -> use bounding box
%
% See also: bwlabeln, bwlabel

dim = ndims(label);
isize = size(label);
labs = imlabel(label);
nextlabel = length(labs) + 1;
seplabel = label;

bb = imlabelboundingboxes(label);

i = 1;
for l = labs
   subim = imextract(label, bb(i,:));
   cc = bwconncomp(subim == l);
   for nl = 2:cc.NumObjects
      
      if dim == 2
         [ix, iy] = ind2sub(size(subim), cc.PixelIdxList{nl});
         ix = ix + bb(i, 1) -1;
         iy = iy + bb(i, 2) -1;
         seplabel(ix, iy) = nextlabel;
      else
         [ix, iy, iz] = ind2sub(size(subim), cc.PixelIdxList{nl});
         ix = ix + bb(i, 1) -1;
         iy = iy + bb(i, 2) -1;
         iz = iz + bb(i, 3) -1;
         seplabel(ix, iy, iz) = nextlabel;
      end

      nextlabel = nextlabel + 1;
   end
   i = i + 1;
end

end
      
      
      