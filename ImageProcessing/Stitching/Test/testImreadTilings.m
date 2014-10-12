%%% compare tiling and image orientation

k = 1; clear imgs
for i = [1,2,3,4, 28, 29, 30, 31, 55, 56, 57, 58]
   fn = tagexpr2string('./Test/Images/hESCells_Tiling_Large/11_150_p<tile,6>t00000001z001c01.tif', 'tile', i)
   imgs{k} = imread(fn)';
   %imgs{k} = impqlpermute(imgs{k}, 'pq', 'py');
   k = k + 1;
end
length(imgs)
imgs = reshape(imgs, [4,3]);

if verbose
   figure(1); clf;
   implottiling(imgs)
end



%%


k = 1; clear imgs
for i = [1,2,3,4, 28, 29, 30, 31, 55, 56, 57, 58]
   fn = tagexpr2string('./Test/Images/hESCells_Tiling_Large/11_150_p<tile,6>t00000001z001c01.tif', 'tile', i)
   imgs{k} = imreadBF(fn);
   imgs{k} = impqlpermute(imgs{k}, 'pq', 'py');
   k = k + 1;
end
length(imgs)
imgs = reshape(imgs, [4,3]);

if verbose
   figure(2); clf;
   implottiling(imgs)
end
