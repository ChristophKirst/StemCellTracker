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

imgsep = imglab;

n = max(imglab(:));
if n == 0
   return
end

dim = ndims(imglab);
isize = size(imglab);
bb = imlabelboundingboxes(imglab);

nextlabel = n +1;

for l = 1:n
   if mod(l, 500) == 0
      fprintf('imlabelseparate: %g / %g\n', l, n);
   end
   %bb(i,:)
   
   bbox = bb(l,:);
   
   % check if we have object at all
   if prod(bbox((dim+1):end)) == 0
      continue
   end

   %flat objects with trailling size 1 do not work with matlab !
   if bbox(dim) == bbox(2*dim)
      if bbox(2*dim) < isize(dim)
         bbox(2*dim) = bbox(dim) + 1;
      else
         bbox(dim) = bbox(dim) - 1;
         if bbox(dim) < 1
            bbox(dim) = 1;
         end
      end
   end
 
   subim = imextract(imglab, bbox);
   %size(subim)
   
   if dim == 2
      conn = 4;
   else
      conn = 18;
   end
   
   cc = bwconncomp(subim == l, conn);
   
   for nl = 2:(cc.NumObjects)
      sub = imind2sub(size(subim), cc.PixelIdxList{nl});
      sub = sub + repmat(bbox(1:dim)-1, size(sub,1), 1);
      ind = imsub2ind(isize, sub);
      imgsep(ind) = nextlabel;
      nextlabel = nextlabel + 1;
   end
end

end
      
      
      