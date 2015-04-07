%%

initialize
bfinitialize


verbose = true;

%%
clear all
clear classes
close all
clc
verbose = true;
initializeParallelProcessing;



%% Setup Image Source

% infer the tag expression o he files from the file folder automatically
%texp = ('/Users/fetoc/Desktop/data/dilutions experiments xBMP4 91914/dilution PB31-xBMP4 c5 in RUES2 1100 fates after 48hrs Dox/Sox2 Bra Myc/20140914_ESCEpFgf52.zvi')
%texp='exp<N>p<tile,2>c<ch>.tif'
texp='1_p<P,6>t00000001z001c<C,2>.tif'

texp = '/var/run/media/ckirst/ChristophsData/for CK_in_situ/1_p<P,6>t00000001z001c<C,2>.tif';

is = ImageSourceFiles(texp);
is.printInfo

%% 

is.setReshape('P','UV',[13,18]);
is.setCellFormat('UvC')
is.setRange('C',1, 'U', 1:3, 'V', 1:3)

is.printInfo

% %%
% 
% im=is.data(1);
% figure(4);
% implot(im)

%% Dangerous for large dataset, Restrict the range first!

figure(1);
is.plottiling


%% Alignment

% create 
clc;
isalgn = Alignment(is);
isalgn.printInfo

%var2char({[algn.pairs.from], [algn.pairs.to]})

%% Background intensity for testing overlap quality of tiles

img1 = is.data(2);
nbins = 50;
th = thresholdFirstMin(img1, 'nbins', nbins, 'delta', 1/1000 * numel(img1))
%th=2000
if verbose
   figure(3)
   hist(img1(:), nbins)
end

%% Quality of Overlap between neighbouring tiles 

% parameter: see overlapQuality
isalgn.calculateOverlapQuality('threshold.max', th, 'overlap.max', 120);
[isalgn.pairs.quality];


%% Connected Components based on overlap quality

clc
subalgn = isalgn.connectedComponents('threshold.quality', -eps);
nsubalgn = length(subalgn)

if verbose 
%    for i = 1:nsubalgn
%       subalgn(i);
%    end   
   var2char({subalgn.nodes})
end


%% Align components

% clear subalgn
% subalgn = isalgn.connectedComponents('threshold.quality', -eps);
% nsubalgn = length(subalgn)

for s = 1:nsubalgn
    s
   subalgn(s).align('alignment', 'RMS', 'overlap.max', 120, 'overlap.min', 40, 'shift.max', 140);
   if verbose && s < 20
      subalgn(s).printInfo%
      figure(4+s)
%       imsubplot(2, 1, s);
      %subalgn(s).plot('method','Min')
      %implot(subalgn(s).data)
            subalgn(s).plotAlignedPreview('scale', 0.05)
   end
end


%% detect colonies in the aligned images
clc
colonies = [];

 for s = 1:nsubalgn
 %   for s = 2:2
    s
    
    var2char(subalgn(s).imageShifts)
    var2char(subalgn(s).absoluteShiftsAndSize)
    
    figure(5+s); clf
    %implot(subalgn(s).data)
    
    
    %rois = findROIsByOpening(subalgn(s).data, 'threshold', th, 'output','ROIs', 'plot', true, 'strel', 50);
    rois = detectROIsByPeakVolume(subalgn(s), 'radius', 100, 'dilate', 50, 'center', true, 'hmax', th, 'plot', true);
    
    for r = 1:length(rois)
    r
        colonies = [colonies, Colony(subalgn(s), rois{r})]; %#ok<AGROW>
    end
 end
 
 save colonies '-v7.3'

%%

clear all
load colonies
%%
for i=61:69
    
    figure(i)
a=colonies(i).data;
imshow(imadjust(uint16(a)));

end
%%
goodCol500=[4,7,13,18,19,25,26,29,32,35,37,39,40,41,42,44,45,46,47,48,49,50,55,58];
goodCol1000=[2,24,52,59,60,61,67];
save goodCol500 '-v7.3';
save goodCol1000 '-v7.3';

imwrite(imadjust(uint16(a)),'forlearning4.tif','tif')

%%
clear all 
load colonies
col=4
% a=colonies(col).data;
% imshow(imadjust(uint16(a)));

is.setTagRange('ch',{'04'});
is.print
a=colonies(col).data;
figure(2)
imshow(imadjust(uint16(a)));