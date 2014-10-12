function dataFrmt = imfrmtReshapeFormat(rawFrmt, reshapeFrom, reshapeTo)
%
% dataFrmt = imfrmtReshapeFormat(rawFrmt, reshapeFrom, reshapeTo)
%
% description: 
%     returns data format obtained by replacing the reshape formats in the raw formats
%
% input:
%     rawFrmt     format of raw data
%     reshapeFrom use these dimensions for reshaping
%     reshapeTo   use these dimensions to reshape to
%
% output:
%     dataFrmt    data format
%
% note:
%    inversion of coordinate axes is ignored

if ~iscell(reshapeFrom)
   reshapeFrom = {reshapeFrom};
end
if ~iscell(reshapeTo)
   reshapeTo = {reshapeTo};
end

n = length(reshapeFrom);

if n ~= length(reshapeTo)
   error('imfrmtReshape: in consistent input!');
end

if n == 0
   dataFrmt = rawFrmt;
   return
end

% reshape sequentially

dataFrmt  = rawFrmt;

for r = 1:n
   id = find(dataFrmt == reshapeFrom{r}(1), 1);
   if ~isempty(id)
      dataFrmt(ismember(dataFrmt, reshapeFrom{r})) = [];
      dataFrmt = [dataFrmt(1:id-1), reshapeTo{r}, dataFrmt(id:end)];
   end
end

end