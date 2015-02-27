%% ijplot5d large data

initialize
ijinitialize

%%

img = syntheticLabeledImage([200, 200, 10], [15, 15, 4], 60);
imgc = imcolorize(img);

figure(1);   clf;
implot(img)

%min(imgc(:))
%max(imgc(:))

%%
tic
imgmovie = repmat(imgc, 1, 1, 1, 1, 10);
ijplot5d(imgmovie, 'class', 'uint8')
toc


%% via files
 
tic
imgmovie = repmat(imgc, 1, 1, 1, 1, 5);
ijplot5d(imgmovie, 'class', 'uint8', 'write', true)
toc
