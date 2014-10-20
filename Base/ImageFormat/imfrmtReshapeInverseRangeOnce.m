function rawRange = imfrmtReshapeInverseRangeOnce(dataSize, dataFrmt, rawFrmt, reshapeFrom, reshapeTo, reshapeSize, range)
%
% rawRange = imfrmtReshapeOnceInverseRange(dataSize, dataFrmt, rawFrmt, reshapeFrom, reshapeTo, reshapeSize, range)
%
% description: 
%     converts a range for data to the range needed to get the data from the raw format
%
% input:
%     lots of inputs...
%
% output:
%     rawRange       the range to read the raw data
%     rawReshape     the reshape size to use form raw to final data reshaping
% 
% See also: imfrmtReshapeCellData


% reshapeFrom must be 1d
if length(reshapeFrom) ~= 1
   error('imfrmtReshapeInverseRangeOnce: raw reshape dimensio not 1d!');
end

% position of reshaped format in data format
posDataReshapeTo = imfrmtPosition(dataFrmt, reshapeTo);
if length(reshapeTo) ~= length(posDataReshapeTo)
   error('imfrmtReshapeInverseRange: inconsistent inputs!');
end

% position of reshape format in raw format
posRawReshapeFrom = imfrmtPosition(rawFrmt, reshapeFrom);
if length(reshapeFrom) ~= length(posRawReshapeFrom)
   error('imfrmtReshapeInverseRange: inconsistent inputs!');
end

% non reshaped coords
nonReshapeDataFrmt = dataFrmt;
nonReshapeDataFrmt(posDataReshapeTo) = [];

nonReshapeRawFrmt = rawFrmt;
nonReshapeRawFrmt(posRawReshapeFrom) = [];


if sum(ismember(lower(nonReshapeDataFrmt), lower(nonReshapeRawFrmt))) ~= length(nonReshapeDataFrmt)
   error('imfrmtReshapeInverseRange: inconsistent inputs!');
end

% non reshaped coords transform directly
% dataFrmt
% nonReshapeRawFrmt
% dataSize
% range

rawRange = imfrmtReformatRange(dataSize, dataFrmt, nonReshapeRawFrmt, range);

% if no reshape ranges specifed we are done 
if ~any(ismember(lower(fieldnames(range)), num2cell(lower(reshapeTo))))
   return
end

% transfrom reshaped specifications, mapping is trivial in index coords

reshapeDataSize = dataSize(posDataReshapeTo);
reshapeDataFrmt = dataFrmt(posDataReshapeTo);

reshapeToSize = imfrmtReformatSize(reshapeDataSize, reshapeDataFrmt, reshapeTo);

% index before reshaping
rawIdx = imfrmtRangeToIndex(reshapeToSize, reshapeTo, range);

% reshaping is trivial

% assign to reshape from frmt
if isemptystruct(rawRange)
   rawRange = struct(reshapeFrom, rawIdx(:)');
else
   rawRange.(reshapeFrom) = rawIdx(:)';
end

end








   