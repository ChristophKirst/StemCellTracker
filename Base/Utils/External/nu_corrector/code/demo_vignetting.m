function demo_vignetting()
% demo for vignetting correction

%%
% setting
addpath(genpath('.'));

% read and process the image

imgr = mat2gray(img);

tic
[im_vign_corrected,im_vignetting]=vignCorrection_nonPara(round(255 * imgr), );%you may need to input the second parameter and adapt it to the vignetting in the input image
toc

% show results
figure; 
imsubplot(3,1,1); implot(imgr); title('Given Image');
imsubplot(3,1,2); implot(mat2gray(im_vign_corrected)); title('Vignetting-Corrected Image')
imsubplot(3,1,3); implot(mat2gray(im_vignetting)); title('Estimated Vignetting')
