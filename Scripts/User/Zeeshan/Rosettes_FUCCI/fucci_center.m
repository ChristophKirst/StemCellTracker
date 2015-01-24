function rc = fucci_center(img, varargin)

if nargin < 2
   verbose = false;
else
   verbose = varargin{1};
end

%% Find Center

% smooth and find minimum

imgs = mat2gray(filterGaussian(img, 200, [], [], 'replicate'));
%imgs = iminvert(imgs);
imgc = imextendedmin(imgs, 0.01);
imgcl = bwlabeln(imgc);
imgcl = imclearborder(imgcl);


% detect the center
imgss = iminvert(imgs);
%imgss = filterDisk(imgs, 300, 5, -1, -1, 1.0);
%imgs = imdilate(imgs, strel('disk', 400));

imgc = imextendedmax(mat2gray(imgss), 0.05);

imgcl = bwlabeln(imgc);

cstats = imstatistics(imgcl, {'Centroid', 'MaxIntensity'},  imgss);
rcenter = [cstats.Centroid];

rcenter = rcenter(:, and(rcenter(1,:) > 750, rcenter(2,:) < 1250));
rcenter = rcenter(:, and(rcenter(2,:) < 800, rcenter(2,:) > 200));

if isempty(rcenter)
   error('cannot detect rcenter!')
elseif length(rcenter(1, :)) > 1
   rcenter = rcenter(:, 1);
end

rc = rcenter;

if verbose 
   figure(31); clf
   implottiling({img, imgs; imgss, imoverlay(img, imgcl)})
   hold on
   plot(rcenter(1), rcenter(2), 'b*', 'MarkerSize', 25, 'LineWidth', 2)
end
