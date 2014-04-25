 
%% Test Background Removal

lab = syntheticLabeledImage([100, 100], 20);

mm = lab > 0;
back = mat2gray(fspecial2('sphere', [100, 100]));
img = mm + back;
max(back(:))

figure(1)
implottiling({lab, back, img})



%% remove trend in illumination via morphological opening
background = imopen(img,strel('disk',15));
imgpre = img - background;
imgpre = imclip(imgpre,0);

figure(2)
implottiling({img, imgpre})


