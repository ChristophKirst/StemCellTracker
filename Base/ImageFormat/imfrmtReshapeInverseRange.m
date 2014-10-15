function [rawRange, rawReshape] = imfrmtReshapeInverseRange(dataSize, dataFrmt, rawFrmt, reshapeFrom, reshapeTo, reshapeSize, range)
%
% [rawRange, rawReshape] = imfrmtReshapeInverseRange(dataSize, dataFrmt, rawFrmt, reshapeFrom, reshapeTo, reshapeSize, range)
%
% description: 
%     converts a tag range for data to the range needed to get the data from the raw format
%
% input:
%     lots of inputs...
%
% output:
%     rawRange       the range to read the raw data
%     rawReshape     the reshape size to use form raw to final data reshaping
% 
% note:
%     ranges cannot be inversed in general as spcifications might not be separable !
%     for image data usually the format is expanded from one dim to several 
%     in which case the inversion in range is possible. we produce an error if reshapeFrom is not a single character
% 
% See also: imfrmtReshapeCellData

if ~iscell(reshapeFrom)
   reshapeFrom = {reshapeFrom};
end
if ~iscell(reshapeTo)
   reshapeTo = {reshapeTo};
end
if ~iscell(reshapeSize)
   reshapeSize = {reshapeSize};
end

n = length(reshapeFrom);

if n ~= length(reshapeTo) || n ~= length(reshapeSize)
   error('imfrmtReshapeInverseRange: in consistent input!');
end


if isemptystruct(range)
   rawRange = range;
   rawReshape  = reshapeSize;
   return
end

if n == 0
   rawRange = imfrmtReformatRange(dataSize, dataFrmt, rawFrmt, range);
   rawReshape  = reshapeSize;
   return
end

% determine intermediate formats
tempFrmt  = cell(1,n+1);
tempFrmt{1} = rawFrmt;

for r = 1:n-1  
   tempFrmtOut = tempFrmt{r};
   tempFrmtOut(ismember(tempFrmt{r}, reshapeFrom{r})) = [];
   tempFrmt{r+1} = [tempFrmtOut, reshapeTo{r}];   
end

tempFrmt{n+1} = dataFrmt;

tempSize = cell(1,n+1);
tempSize{n+1} = dataSize;
for r = n:-1:1
   tempSize{r} = imfrmtReshapeInverseSizeOnce(tempSize{r+1}, tempFrmt{r+1}, tempFrmt{r}, reshapeFrom{r}, reshapeTo{r}, reshapeSize{r});
end
% 
% disp Reshape
% var2char(tempFrmt)
% var2char(tempSize)
% dataSize


% inverse reshaping the ranges
rawRange = range;
rawReshape = cell(1,n);

for r = n:-1:1
   if nargout > 1
      si = imfrmtReformatSize(tempSize{r+1}, tempFrmt{r+1}, reshapeTo{r});
      rawReshape{r} = imfrmtRangeSize(si, reshapeTo{r}, rawRange);
   end
   
   rawRange = imfrmtReshapeInverseRangeOnce(tempSize{r+1}, tempFrmt{r+1}, tempFrmt{r}, ...
                                               reshapeFrom{r}, reshapeTo{r}, reshapeSize{r}, rawRange);
end

end








   