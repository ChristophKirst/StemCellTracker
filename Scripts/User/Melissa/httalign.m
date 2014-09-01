function shall = httalign(imgs, verbose)


%% Stich


% algin  upper two
movl = 180;

img1 = imgs{1,1}; img2 = imgs{1,2};

img1r = img1(end-movl+1:end, :); img2r = img2(1:movl, :);

%figure(3);
%implottiling({img1r, img2r})

sh = alignGlobally2ImagesByRMS(img1r, img2r);
sh{2} = sh{2} + [size(img1, 1)-movl, 0];

shall = sh;

% align upper lower

img1 = imgs{1,1}; img2 = imgs{2,1};

img1r = img1(:, 1:movl); img2r = img2(:, end-movl+1:end);

%figure(3);
%implottiling({img1r, img2r})

sh = alignGlobally2ImagesByRMS(img1r, img2r);
sh{2} = sh{2} - [0, size(img1,2) - movl] ;


%figure(2); clf;
%plotAlignedImages({mat2gray(img1), mat2gray(img2)}, sh)

shall{2,1} = sh{2};

% align upper lower right

img1 = imgs{1,2}; img2 = imgs{2,2};

img1r = img1(:, 1:movl); img2r = img2(:, end-movl+1:end);

%figure(3);
%implottiling({img1r, img2r})

sh = alignGlobally2ImagesByRMS(img1r, img2r);
sh{2} = sh{2} - [0, size(img1,2) - movl] ;

%figure(2); clf;
%plotAlignedImages({mat2gray(img1), mat2gray(img2)}, sh)

shall{2,2} = sh{2} + shall{1,2};

% final image

%%

if verbose
   figure(2); clf;
   plotAlignedImages(cellfun(@mat2gray, imgs, 'UniformOutput', false), shall)
end

