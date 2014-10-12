function data = imfrmtReshape(rawData, rawFrmt, dataFrmt, reshapeFrom, reshapeTo, reshapeSize)
%
% data = imfrmtReshape(rawData, rawFrmt, dataFrmt, reshapeFrom, reshapeTo, reshapeSize)
%
% description: 
%    reshapes dims reshapeFrom of data of foramt inFrmt to format reshapeTo 
%    using the reshape size reshapeSize and assuming the full output format outFrmt
%
% input:
%     rawData     raw input data or cell
%     rawFrmt     format of raw data
%     dataFrmt    format of output data
%     reshapeFrom use these dimensions for reshaping
%     reshapeTo   use these dimensions to reshape to
%     reshapeSize use this size for reshaping
%
% output:
%     data        reshape data

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
   error('imfrmtReshape: in consistent input!');
end

if n == 0
   data = imfrmtReformat(rawData, rawFrmt, dataFrmt);
   return
end

% reshape sequentially

activeInFrmt  = rawFrmt;
activeOutFrmt = rawFrmt;
data = rawData;

for r = 1:n-1
   activeOutFrmt(ismember(activeOutFrmt, reshapeFrom{r})) = [];
   activeOutFrmt = [activeOutFrmt, reshapeTo{r}]; %#ok<AGROW>

   data = imfrmtReshapeOnce(data, activeInFrmt, activeOutFrmt, ...
                               reshapeFrom{r}, reshapeTo{r}, reshapeSize{r});
   
   activeInFrmt = activeOutFrmt;
end

data = imfrmtReshapeOnce(data, activeInFrmt, dataFrmt,...
                               reshapeFrom{end}, reshapeTo{end}, reshapeSize{end});

end