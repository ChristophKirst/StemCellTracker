function range = imfrmtReshapeRange(rawSize, rawFrmt, dataFrmt, reshapeFrom, reshapeTo, reshapeSize, range)
%
% data = imfrmtReshape(rawData, rawFrmt, dataFrmt, reshapeFrom, reshapeTo, reshapeSize)
%
% description: 
%    reshapes a range
%
% input:
%     rawSize     raw size of input data
%     rawFrmt     format of raw data
%     dataFrmt    format of output data
%     reshapeFrom use these dimensions for reshaping
%     reshapeTo   use these dimensions to reshape to
%     reshapeSize use this size for reshaping
%     range       the range to reshape
%
% output:
%     outRange    reshaped range

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
   range = imfrmtReformatRange(rawSize, rawFrmt, dataFrmt, range);
   return
end

% reshape sequentially

activeInFrmt  = rawFrmt;
activeOutFrmt = rawFrmt;

for r = 1:n-1
   activeOutFrmt(ismember(activeOutFrmt, reshapeFrom{r})) = [];
   activeOutFrmt = [activeOutFrmt, reshapeTo{r}]; %#ok<AGROW>

   range = imfrmtReshapeRangeOnce(rawSize, activeInFrmt, activeOutFrmt, ...
                                  reshapeFrom{r}, reshapeTo{r}, reshapeSize{r}, range);
                               
   rawSize = imfrmtReshapeSizeOnce(rawSize, activeInFrmt, activeOutFrmt, ...
                                  reshapeFrom{r}, reshapeTo{r}, reshapeSize{r});
   
   activeInFrmt = activeOutFrmt;
end

range = imfrmtReshapeRangeOnce(rawSize, activeInFrmt, dataFrmt,...
                               reshapeFrom{end}, reshapeTo{end}, reshapeSize{end}, range);

end