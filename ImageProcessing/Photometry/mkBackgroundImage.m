function [minIm, meanIm]=mkBackgroundImage(direc,fileKeyword,maxIm,filterrad)


if ~exist('filterrad','var')
    filterrad=200;
end

[~, ImFiles]=folderFilesFromKeyword(direc,fileKeyword);

if isempty([ImFiles.name])
    disp('folderFilesFromKeyword failed, using dir...');
    ImFiles=dir([direc filesep '*' fileKeyword '*']);
end

if isempty([ImFiles.name])
    minIm=[]; meanIm=[];
    return;
end

q=1;
nIms=length(ImFiles);
ImRange=randperm(nIms-1);
for jj=1:length(ImRange)
    ii=ImRange(jj);
    if exist('maxIm','var') && q > maxIm
        break;
    end
    imnm=ImFiles(ii).name;
    imNow=im2double(imread([direc filesep imnm]));
    %if max(max(imNow)) > 300
        if q==1
            minIm=imNow;
            meanIm=imNow;
        else
            minIm=min(minIm,imNow);
            meanIm=((q-1)*meanIm+imNow)/q;
        end
        q=q+1;
%      else
%          continue;
%      end
end

 gfilt=fspecial('gaussian',filterrad,filterrad/5);
 minIm=imfilter(minIm,gfilt,'symmetric');
 meanIm=imfilter(meanIm,gfilt,'symmetric');