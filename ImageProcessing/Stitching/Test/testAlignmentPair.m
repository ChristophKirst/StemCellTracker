%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Image Alignment Class %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
clear classes
close all

initialize


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AlignmentPair calss

%% construction and routines
ap = AlignmentPair(1,2,'lr')
ap.dim


c = ap.toCell()


%%
ap = AlignmentPair(1,2,'du')

c = ap.toCell()

ap2 = AlignmentPair(c)

%%

ap = AlignmentPair(1,2,'bt')

c = ap.toCell()


%% alignImagePair

%% Load Data
img = loadTestImage();
img = mat2gray(img);
size(img)

img1 = img(1:230,1:end-15);
img2 = img(200:end,20:end);

figure(2);
implottiling({img1;img2}, 'link', false);


%%
p = AlignmentPair(img1, img2, 'lr')

%%
p = alignImagePair(p, 'shift.max', 50)

figure(3); clf
p.plot()

%%

tic
p.align('shift.max', 50, 'overlap.max', 50, 'overlap.min', 20)
toc

figure(3); clf
p.plot()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Alignment classes


%% constructor

igird = {1,2; 3,4; 5,6}';

a = Alignment(igird)

[a.pairs.from]
[a.pairs.to]
[a.pairs.orientation]

%% 


