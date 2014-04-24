function imgsep = imlabelseparate(imglab)
%
% imgsep = imlabelseparate(imglab)
%
% description:
%     relabels the image by separating non touching regions with same imglab
%     into multiple imglab.
%
% input:
%     imglab       labeled image
% 
% output:
%     imgsep       relabled image, no separate regoins with same label
%
% note:
%     bwlabeln(imglab > 0) fails if two labels touch
%     bwconncomp scales badly with image size even for small objects -> use bounding box
%
% See also: bwlabeln, bwlabel

dim = ndims(imglab);
%isize = size(imglab);
labs = imlabel(imglab);
nextlabel = length(labs) + 1;
imgsep = imglab;

bb = imlabelboundingboxes(imglab);

n = length(labs);
i = 1;
for l = labs
   if mod(i, 50) == 0
      fprintf('imlabelseparate: %g / %g\n', i, n);
   end
   %bb(i,:)
   
   subim = imextract(imglab, bb(i,:));
   %size(subim)
   
   if dim == 2
      conn = 4;
   else
      conn = 18;
   end
   
   cc = bwconncomp(subim == l, conn);
   
   for nl = 2:(cc.NumObjects)
      

      %dim
      %size(subim)
      %cc.PixelIdxList{nl}
      %bb(i,:)

      if dim == 2
         [ix, iy] = ind2sub(size(subim), cc.PixelIdxList{nl});
         ix = ix + bb(i, 1) -1;
         iy = iy + bb(i, 2) -1;
         
         %[ix,iy]
         
         imgsep(ix, iy) = nextlabel;
      else
         [ix, iy, iz] = ind2sub(size(subim), cc.PixelIdxList{nl});
         ix = ix + bb(i, 1) -1;
         iy = iy + bb(i, 2) -1;
         iz = iz + bb(i, 3) -1;
         
         %[ix,iy,iz]
         %length(cc.PixelIdxList{nl}) 
         %{min(ix), max(ix)}
         %{min(iy), max(iy)}
         %{min(iz), max(iz)} 
         %size(imgsep)
         %nextlabel
         
         for k  = 1:length(ix) % strange matlab bug???
            %[k ix(k), iy(k), iz(k)]
            imgsep(ix(k), iy(k), iz(k)) = nextlabel;
         end
      end

      nextlabel = nextlabel + 1;
   end
   i = i + 1;
   
   drawnow
end

end
      
      
      