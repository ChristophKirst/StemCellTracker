function [dataFrmt, cellFrmt] = imfrmtReshapeCellDataFormat(rawDataFrmt, rawCellFrmt, reshapeFrom, reshapeTo)
%
% [dataFrmt, cellFrmt] = imfrmtReshapeCellDataFormat(rawDataFrmt, rawCellFormat, reshapeFrom, reshapeTo)
%
% description: 
%     returns data and cell format obtained by replacing the reshape formats in the raw formats
%
% input:
%     raw*Frmt    format of raw data / cell
%     reshapeFrom use these dimensions for reshaping
%     reshapeTo   use these dimensions to reshape to
%
% output:
%     *Frmt       reshaped data/cell format


dataFrmt = imfrmtReshapeFormat(rawDataFrmt, reshapeFrom, reshapeTo);
cellFrmt = imfrmtReshapeFormat(rawCellFrmt, reshapeFrom, reshapeTo);

end