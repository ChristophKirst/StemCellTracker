%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Image Formats    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all
initialize


%% Load data

img =loadTestImage;
img = img(:,1:300);

figure(1)
implot(img)

%% imfrmtFormat

clc
imfrmtFormat(img)


%% 

clc

imfrmtFormat(cell(1,3))
imfrmtFormat(cell(3,1))


%% imfrmtPosition

imfrmtPosition('XYZC', 'yxce')


%% imfrmtSpatialDims

imfrmtSpatialDims('XYZC')
imfrmtSpatialDims('XCY')

%% imfrmtAssignment

asgn = imfrmtAssignment([100,1, 200], 'XZY', 'X', 1:2, 'y', 1:5);
var2char(asgn)

%% imfrmtDataSubset

imgs  = imfrmtDataSubset(img, 'XY', 'X', 1:100, 'y', 1:100);
imgs2 = imfrmtDataSubset(img, 'XY', 'X', 1:100, 'Y', 1:100);

figure(1); clf 
implottiling({img; imgs; imgs2}, 'link', false)

 
%% imfrmtReformat  

imgr = imfrmtReformat(img, 'XY', 'yX');

figure(1); clf 
implottiling({img; imgr}, 'link', false)


%%
% data loss is accepted
imgr = imfrmtReformat(img, 'XY', 'X')
img(:,1) - imgr


%% imfrmtReformatSize

imfrmtReformatSize([10,20,30], [], 'XCZy')

%%

%data loss is accepted
imfrmtReformatSize([1,2,3,4,5], 'XYZCT', 'zX')


%% imfrmtPermute

clc
imgYX = imfrmtPermute(img, 'XY', 'YX');
imgXy = imfrmtPermute(img, 'XY', 'Yx');

figure(2); clf
implottiling({img; imgYX; imgXy}, 'link', false)


%%

% extension adds extra dim
imgXZY = imfrmtPermute(img, 'XY', 'XZY');
size(imgXZY)


%% 

% loosing data produces error
imgError = imfrmtPermute(img, 'XY', 'X');
size(imgXZY)


%% imfrmtPermuteSize

clc
imfrmtPermuteSize(size(img), 'XY', 'XY')
imfrmtPermuteSize(size(img), 'XY', 'yx')
imfrmtPermuteSize(size(img), 'XY', 'yxc')

%%

% losing data produces error
imfrmtPermuteSize(img, 'XY', 'X');



%% Cells


c = num2cell(reshape(1:20, 4, 5))
imfrmtReformat(c, 'UV', 'uV')


%% Ranges


tgr.X = 1:10;
tgr.Y = 1:5;


imfrmtReformatRange(tgr, 'XY', [100,200], 'y')


%% Combined Cell + Data arrays
clc
for i = 1:4
   for j = 1:7
      c{i,j} = rand(3,2);
   end
end

cr = imfrmtReformatCell(c, 'UV', 'XY', 'uV', 'XY');
size(cr)

a = 1; b = 2;
cr{a,b}
c{end,b}
c{a,b}


%%
clc

cr = imfrmtReformatCell(c, 'UV', 'XY', 'UX', 'VY');
size(cr)
size(cr{1})

a = 1; b = 2;
cr{a,b}
c{a,b}



%% imfrmtReformatIndex

img = rand(10,30);

id = [5,8,14,20]';

sid = imind2sub(size(img), id)


idrr = imfrmtReformatIndex(id, size(img), 'XY', 'yX');
[sidr, isr] = imfrmtReformatSubIndex(sid, size(img), 'XY', 'yX');
isr

imgr = imfrmtReformat(img, 'XY', 'yX');

idr = imsub2ind(isr, sidr)

img(id)'
imgr(idr)'
imgr(idrr)'




%% imfrmtRangeSize

imfrmtRangeSize([10,3,8,6], 'XYCZ', setParameter('X', [1:5], 'C', 1:3))





