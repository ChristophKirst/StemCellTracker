function [im_corrected,im_vignetting]=vignCorrection_nonPara(im_given,varargin)
% % 
% % vignCorrection_nonPara can esitmate and remove vignetting effect in gray
% % or color images. Vignetting image is represented nonparamatrically as
% % shown in the reference paper.
% % 
% % Examples:
% % vignCorrection_nonPara(im_given)
% % vignCorrection_nonPara(im_given,labmda)
% % vignCorrection_nonPara(im_given,labmda,itrNum)
% % vignCorrection_nonPara(im_given,labmda,itrNum,alpha)
% % 
% % Input Arguments:
% % im_given: Matrix of image data, obtained with imread or other ways.
% % labmda:   Parameter adjusting the smoothness of the resulting vignetting
% %           image. Default is 0.3. See Eq. (16)
% % itrNum:   Number of iterations for the IRLS. Default is 4.
% % alpha:    Exponent value of gradient distribution, in Eq. (14)  
% % 
% % Ouput Arguments:
% % im_corrected:       image after vignetting correction
% % im_vignetting:      image of the estimated vignetting
% % 
% % Reference:
% % Y. Zheng, J. Yu, S. B. Kang, S. Lin, and C. Kambhamettu, �Single-image vignetting
% % correction using radial gradient symmetry,� in IEEE Conference on
% % Computer Vision and Pattern Recognition, June 2008, pp. 1�8.
% % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % Re-implemented by Yuanjie Zheng on April 16, 2010 @PICSL@UPenn,
% % based on the original work for the conference paper.
% % Contact: zheng.vision@gmail.com
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% check input
% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if nargin==0
    error('You need to provide image data for vignCorrection_nonPara.m. Try "help vignCorrection_nonPara" to get more information.');
end

% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% set up parameters
% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% variable paras
[lambda,itrNum,alpha]=parseInputs(nargin-1,varargin);
% fixed paras
epsilon=0.0001;% perturbation on B
shift=1; % shift to aoivd computation of log(0)

% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% downsample to increase speed 
% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dsfact=0.25;
im_given_sampled=imresize(im_given,dsfact);

% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% convert color to gray
% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
numsz=length(size(im_given_sampled));
if numsz==3
    DIM=3;
    im_gray=rgb2gray(uint8(im_given_sampled));
else
    DIM=1;
    im_gray=im_given_sampled;
end
im_data=double(im_gray);

sz=size(im_data);
numPixels=sz(1)*sz(2);


% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% preparitions
% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
disp('Preparing ...')
Z_shift=log(im_data+shift+eps); %%%%%%%%%%important
% Z_orn=log(im_data+eps); 

% initial W to 1 for each pixel
vector_W=ones(numPixels,1);
W=W_vec2sparse(vector_W); 

% compute Lvi
Lvxny=getLv(sz,1)+getLv(sz,2);

% get C
C=getC(sz);
numR=size(C,2);

% compute Lvii
Lxxnyy=getLxx([numR],1);

%  identiy matrix
myI=speye(numR,numR);

Gamma1=epsilon*myI;
Gamma2=lambda*(2*numPixels/numR)*Lxxnyy;

% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Iterative reweighted least squares
% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
colNum=itrNum+1;
ims{1}=im_given; tts{1}='Given Image';
ims{colNum+1}=uint8(ones(sz(1),sz(2))*255); tts{colNum+1}='Intial Vignetting';
ims{2*colNum+1}=uint8(reshape(vector_W,sz)*255); tts{2*colNum+1}='Intial Weight';
for i=1:itrNum
    disp(['IRLS: iteration ' num2str(i) ' ... ...'])
    
    %solve B
    right=W*Lvxny*Z_shift(:);
    A=W*Lvxny;
    B_r=((C'*A'*A*C+Gamma1'*Gamma1+Gamma2'*Gamma2))\(C'*A'*right);
    B=C*B_r;
    
    % update W
    S1=abs(Lvxny*B(:)-Lvxny*Z_shift(:)); S2=alpha*(S1).^(alpha-1);
    vector_W=exp(-S1).*(1-exp(-S2));    
    W=W_vec2sparse(vector_W);
    
    image_W=uint8(reshape(vector_W,sz)*255);
    
    % 
    B=B-mean(B(:));%this can make the mean value of the corrected image equals to the one of the original image, this also makes the estimated bias's values around 1
    
    %get image for showing estimated vignetting
    b=exp(B-(max(B(:))));% such that b is within [0 1], for display
    image_b=reshape(b,sz);
    image_b=uint8(image_b*255);
    
    % get corrected image
    X=log(double(im_given_sampled)+shift)-repmat(reshape(B,sz(1),sz(2)),[1 1 DIM]);
    x=exp(X)-shift;
    x=uint8(x);
    
    % put all images in a cell array in order to display them together
    ims{1+i}=x; tts{1+i}=['Corrected: Iter ' num2str(i)];
    ims{colNum+1+i}=image_b; tts{colNum+1+i}=['Vignetting: Iter ' num2str(i)];
    ims{2*colNum+1+i}=image_W; tts{2*colNum+1+i}=['Weight: Iter ' num2str(i)];
end


% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% upsample
% % % % % % % % % % % % % % % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% bias
im_vignetting=imresize(image_b,[size(im_given,1) size(im_given,2)],'bicubic');

% bias corrected image
im_bias_4compute=reshape(exp(B),sz)*255;
im_bias_4compute=imresize(im_bias_4compute,[size(im_given,1) size(im_given,2)],'bicubic');

data_bias=double(im_bias_4compute)/255;
X=log(double(im_given)+shift)-repmat(log(data_bias),[1 1 DIM]);
x=exp(X)-shift;
im_corrected=uint8(x);

function [lambda,itrNum,alpha]=parseInputs(numVar,vars)
% default
lambda=0.3;% smoothness
itrNum=4;
alpha=0.6;%diff 2 weight %smaller alpha means more contingent or sharper

% read input
if numVar>=1
    lambda=vars{1};
end
if numVar>=2
    itrNum=vars{2};
end
if numVar>=3
    alpha=vars{3};
end

