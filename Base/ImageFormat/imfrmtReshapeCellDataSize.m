function [outDataSize, outCellSize] = imfrmtReshapeCellDataSize(inDataSize, inCellSize, inDataFrmt, inCellFrmt, outDataFrmt, outCellFrmt, reshapeFrom, reshapeTo, reshapeSize)
%
% outData = imfrmtReshapeCellData(inData, inCellFrmt, inDataFrmt, outCellFrmt, outDataFromat, reshapeFrom, reshapeTo, reshapeSize)
%
% description: 
%     reshapes size as the result from imfrtReshapeCellData 
%
% input:
%     inDataSize  size of data
%     inCellSize  size of cell
%     in*Frmt     format of input image
%     out*Frmt    format of output image 
%     reshapeFrom use these dimensions for reshaping
%     reshapeTo   use these dimensions to reshape to
%     reshapeSize use this size for reshaping
%
% output:
%     outDataSize  size of data
%     outCellSize  size of cell

% convert to full data array
fullSize = [inDataSize, inCellSize];
fullInFrmt = [inDataFrmt, inCellFrmt];
fullOutFrmt = [outDataFrmt, outCellFrmt];

% reshape
fullSize =imfrmtReshapeSize(fullSize, fullInFrmt, fullOutFrmt, reshapeFrom, reshapeTo, reshapeSize);

% split
outDataSize = fullSize(1:length(outDataFrmt));
outCellSize = fullSize((length(outDataFrmt)+1):end);

end