function outSize = imfrmtReshapeSize(inSize, inFrmt, outFrmt, reshapeFrom, reshapeTo, reshapeSize)
%
% outDataSize = imfrmtReshapeSize(inDataSize, inCellSize, inFrmt, outFrmt, reshapeFrom, reshapeTo, reshapeSize)
%
% description: 
%    transforms sizes as imfrmtReshape would do for data
%
% input:
%     inSize      input 
%     infrmt      format of input image
%     outfrmt     format of output image 
%     reshapeFrom use these dimensions for reshaping
%     reshapeTo   use these dimensions to reshape to
%     reshapeSize use this size for reshaping
%
% output:
%     outSize     reshaped size
%
% note:
%    inversion of coordinate axes is ignored

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
   outSize = imfrmtReformatSize(inSize, inFrmt, outFrmt);
   return
end

% reshape sequentially

activeInFrmt  = inFrmt;
activeOutFrmt = inFrmt;
outSize = inSize;

for r = 1:n-1
   activeOutFrmt(ismember(activeInFrmt, reshapeFrom{r})) = [];
   activeOutFrmt = [activeOutFrmt, reshapeTo{r}]; %#ok<AGROW>
   
   
   outSize = imfrmtReshapeSizeOnce(outSize, activeInFrmt, activeOutFrmt, ...
                               reshapeFrom{r}, reshapeTo{r}, reshapeSize{r});
   
   activeInFrmt = activeOutFrmt;
end

outSize = imfrmtReshapeSizeOnce(outSize, activeInFrmt, outFrmt,...
                               reshapeFrom{end}, reshapeTo{end}, reshapeSize{end});

end