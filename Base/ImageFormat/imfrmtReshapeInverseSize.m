function rawSize = imfrmtReshapeInverseSize(dataSize, dataFrmt, rawFrmt, reshapeFrom, reshapeTo, reshapeSize)
%
% outSize = imfrmtReshapeInverseSize(dataSize, dataFrmt, inFrmt, reshapeFrom, reshapeTo, reshapeSize)
%
% description: 
%    transforms sizes of data as imfrmtReshape would do for data
%
% input:
%     dataSize    input data size
%     dataFrmt    format of input 
%     rawFrmt     raw format of raw data
%     reshapeFrom use these dimensions for reshaping
%     reshapeTo   use these dimensions to reshape to
%     reshapeSize use this size for reshaping
%
% output:
%     rawSize     inversely reshaped size
%
% note:
%     as ranges cannot be inverted as long as the reshaping dim is not 1d we restrict to 1d reshapeFrom formats


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
   error('imfrmtReshapeInverseSize: in consistent input!');
end

if n == 0
   rawSize = imfrmtReformatSize(dataSize, dataFrmt, rawFrmt);
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

% inverse reshaping the size
rawSize = dataSize;
for r = n:-1:1
   rawSize = imfrmtReshapeInverseSizeOnce(rawSize, tempFrmt{r+1}, tempFrmt{r}, reshapeFrom{r}, reshapeTo{r}, reshapeSize{r});
end

end
