%% ijplot5d large data

ijinitialize

%%

img = syntheticLabeledImage([512, 512, 20], [50, 50, 4], 60);
imgc = imcolorize(img);

%min(imgc(:))
%max(imgc(:))

tic
imgmovie = repmat(imgc, 1, 1, 1, 1, 10);
ijplot5d(imgmovie, 'class', 'uint8')
toc



%% via files

tic
imgmovie = repmat(imgc, 1, 1, 1, 1, 10);
ijplot5d(imgmovie, 'class', 'uint8', 'write', true)
toc
