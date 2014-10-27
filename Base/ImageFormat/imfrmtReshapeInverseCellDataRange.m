function [rawRange, rawReshape] = imfrmtReshapeInverseCellDataRange(dataSize, dataFrmt, rawDataFrmt, cellSize, cellFrmt, rawCellFrmt, reshapeFrom, reshapeTo, reshapeSize, range)
%
% [rawRange, rawReshape] = imfrmtReshapeInverseCellDataRange(dataSize, dataFrmt, rawFrmt, reshapeFrom, reshapeTo, reshapeSize, range)
%
% description: 
%     converts a range for data to the range needed to get the data from the raw format
%
% input:
%     lots of inputs...
%
% output:
%     rawRange       the coord range to read the raw data
%     rawReshape     the reshape size to use form raw to final data reshaping
% 
% note:
%     ranges cannot be inversed in general as specifications might not be separable !
%     for image data usually the format is expanded from one dim to several 
%     in which case the inversion in range is possible. we produce an error if reshapeFrom is not a single character
% 
% See also: imfrmtReshapeCellData

% combine the formats and sizes

fullSize = [dataSize, cellSize];
fullFrmt = [dataFrmt, cellFrmt];
fullRawFrmt = [rawDataFrmt, rawCellFrmt];

[rawRange, rawReshape] = imfrmtReshapeInverseRange(fullSize, fullFrmt, fullRawFrmt, reshapeFrom, reshapeTo, reshapeSize, range);

end








   