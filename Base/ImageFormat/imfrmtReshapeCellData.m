function outData = imfrmtReshapeCellData(inData, inDataFrmt, inCellFrmt, outDataFrmt, outCellFrmt, reshapeFrom, reshapeTo, reshapeSize)
%
% outData = imfrmtReshapeCellData(inData, inDataFrmt,  inCellFrmt, outDataFromat, outCellFrmt, reshapeFrom, reshapeTo, reshapeSize)
%
% description: 
%     reshapescellas well as 
%
% input:
%     data        input data or cell
%     in*Frmt     format of input image
%     out*Frmt    format of output image 
%     reshapeFrom use these dimensions for reshaping
%     reshapeTo   use these dimensions to reshape to
%     reshapeSize use this size for reshaping
%
% output:
%     outData     reshape data        reformatted image
%
% note:
%    inversion of coordinate axes is ignored


% convert to full data array
outData = imfrmtCellDataToData(inData, inDataFrmt, inCellFrmt);

fullInFrmt = [inDataFrmt, inCellFrmt];
fullOutFrmt = [outDataFrmt, outCellFrmt];

% reshape
outData =imfrmtReshape(outData, fullInFrmt, fullOutFrmt, reshapeFrom, reshapeTo, reshapeSize);

% convert to cell array
outData = imfrmtDataToCellData(outData, outDataFrmt, outCellFrmt);

end