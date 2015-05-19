%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test Image for Non/linear alignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initialize

%%


coords = {2:10:100, 2:10:100, 2:4:15};
coords2 = coords;
for i = 1:3
   coords2{i} = coords2{i} +  randi([-1,1], 1, length(coords2{i}));
end

img1 = zeros(110, 110, 20);
img2 = img1;

%% put sphere in each coord

sphere = fspecial3('sphere', [7,7,7]);
%figure(1); clf; 
%implot(sphere)


%%

ncoordsX = length(coords{1});
ncoordsY = length(coords{2});
ncoordsZ = length(coords{3});

for x = 1:ncoordsX
   for y = 1:ncoordsY
      for z = 1:ncoordsZ
         img1 = imreplace(img1, sphere, [coords{1}(x), coords{2}(y), coords{3}(z)]);
         img2 = imreplace(img2, sphere, [coords2{1}(x), coords2{2}(y), coords2{3}(z)], 'chop', true);
      end
   end
end


%%

figure(6); clf
implottiling(img1, 'tiling', [5,4])


figure(7); clf
implottiling(img2, 'tiling', [5,4])




%%

shifts = XXX;

%% get points from segmentation

pts = XXX;




%% Down Sampling etc test out downstairs.




