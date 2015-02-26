%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test GMM for Segmentation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initialize
bfinitialize

%% Load Image

img  = imreadBF('./Test/Images/hESCells_DAPI.tif');

%img = mat2gray(img(300:370, 240:310));
%img = mat2gray(img(300:340, 245:285));

figure(1); clf
implot(img)

%% Generate Data from Image

% total points
n = 5000;
[x,y] = meshgrid(1:size(img,1), 1:size(img,2));


imgn = fix(n*img/total(img));
figure(1); clf
implottiling({img; mat2gray(imgn)})

data = [];
for i = 1:numel(img)
   data = [data, repmat([x(i); y(i)], 1, imgn(i)) + rand(2, imgn(i))];
end


figure(9); clf
scattercloud(data(2,:), data(1,:), 50, 1)


%% fit GMM

X= numpyFromMat(data');
cp = numpyFromMat(5);
cp  = cp.astype('int');

X.tofile('/home/ckirst/Desktop/testNp.dat')

%%
gmm = py.sklearn.mixture.GMM(pyargs('n_components', cp, 'covariance_type', 'full'))
gmm.fit(X)

%%
clc
%imgf = filterBM(img);

%%

%imglab = segmentBySLIC(img, 'superpixel', 45, 'compactness', 10);
imglab = segmentBySLIC(img, 'superpixel', 2045, 'compactness', 1);

%%
imglab = watershed(1-mat2gray(img));

%%
figure(1); clf


for i = 1:500;
   imgmask = mat2gray(img) > 0.001 * i;
   m(i) = total(imgmask);
end
plot(m)


%%

imgmask = mat2gray(img) > 0.001 * 18;
imgmask = imopen(imgmask, strel('disk', 2));
implot(imgmask)


%%

imgf = filterBM(repmat(mat2gray(img),1,1,3), 'sigma', 10);
imgf = imgf(:,:,1);

figure(4); clf;
implottiling({img; mat2gray(imgf)})

%%
clc
imgf2 = filterSphere(imgf, [3,3]);


figure(5); clf;
implottiling({mat2gray(img); mat2gray(imgf2)})

%%

%%

ds= [2, 3,4, 5,6, 7, 9]

for i = 1:length(ds)
   clc
   imgf2 = filterDisk(imgf, ds(i) * [1,1]);
   imgf2 = imgf - imgf2;
   imgfs{i} = mat2gray(imgf2);
end

imge =mat2gray(imgf)+mat2gray(imgfs{end});
%imge = imge - min(imge(:));
imge = imge - 0.5;
imge(imge<0) = 0;

imgse = imextendedmax(imge, 0.1);


% figure(5); clf;
% hist(imge(:), 256);

figure(6); clf;
implottiling({mat2gray(img), imoverlaylabel(mat2gray(imge), imgse), imgfs{:}}', 'tiling', [4,2])


%%

%imglab = watershed(1-mat2gray(imgf(:,:,1)));
imglab = segmentBySLIC(img, 'superpixel', 2045, 'compactness', 5);
imglab = immask(imglab, imgmask);

figure(2); clf
implottiling({img; imoverlaylabel(img, impixelsurface(imglab), false)})









% %%
% 
%  Fit a mixture of Gaussians with EM using five components
% gmm = mixture.GMM(n_components=5, covariance_type='full')
% gmm.fit(X)
% 
% # Fit a Dirichlet process mixture of Gaussians using five components
% dpgmm = mixture.DPGMM(n_components=5, covariance_type='full')
% dpgmm.fit(X)
% 
% color_iter = itertools.cycle(['r', 'g', 'b', 'c', 'm'])
% 
% for i, (clf, title) in enumerate([(gmm, 'GMM'),
%                                   (dpgmm, 'Dirichlet Process GMM')]):
%     splot = plt.subplot(2, 1, 1 + i)
%     Y_ = clf.predict(X)
%     for i, (mean, covar, color) in enumerate(zip(
%             clf.means_, clf._get_covars(), color_iter)):
%         v, w = linalg.eigh(covar)
%         u = w[0] / linalg.norm(w[0])
%         # as the DP will not use every component it has access to
%         # unless it needs it, we shouldn't plot the redundant
%         # components.
%         if not np.any(Y_ == i):
%             continue
%         plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], .8, color=color)
% 
%         # Plot an ellipse to show the Gaussian component
%         angle = np.arctan(u[1] / u[0])
%         angle = 180 * angle / np.pi  # convert to degrees
%         ell = mpl.patches.Ellipse(mean, v[0], v[1], 180 + angle, color=color)
%         ell.set_clip_box(splot.bbox)
%         ell.set_alpha(0.5)
%         splot.add_artist(ell)
% 
%     plt.xlim(-10, 10)
%     plt.ylim(-3, 6)
%     plt.xticks(())
%     plt.yticks(())
%     plt.title(title)
% 
% plt.show()
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 






%%