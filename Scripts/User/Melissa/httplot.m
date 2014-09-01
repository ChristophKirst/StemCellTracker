function httplot(filename)

%% Plot the Image

imgs = imread_bf(filename, 'series', [1,2, 3, 4], 'channel', 1);

imgs = cellfun(@(x) squeeze(x), imgs, 'UniformOutput', false);

imgs = [imgs(3), imgs(4); imgs(1), imgs(2)];
isizes = cellfun(@size, imgs, 'UniformOutput', false);

sh = httalign(imgs, false)


%% correct illumination etc

for t = 1:4
   imgback = imopen(imgs{t}, strel('disk', 20));
   imgsc{t}  = imgs{t} - imgback;
end
imgsc = reshape(imgsc, 2,2);


%% stitch
img = stitchImagesByMean(imgsc, sh);



%% load other channels and stich

imgch{1} = img;

for c = 2:4
   
   imgs = imread_bf(filename, 'series', [1,2, 3, 4], 'channel', c);
   imgs = cellfun(@(x) squeeze(x), imgs, 'UniformOutput', false);
   imgs = [imgs(3), imgs(4); imgs(1), imgs(2)];

   imgch{c} = stitchImagesByMean(imgs, sh);
end



figure(3); imcolormap('gray');
implottiling(reshape(imgch,[2,2]))


%% f
