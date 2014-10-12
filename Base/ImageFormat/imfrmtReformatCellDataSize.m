function [dataSize, cellSize] = imfrmtReformatCellDataSize(dataSize, cellSize, inDataFrmt, inCellFrmt, outDataFrmt, outCellFrmt)
%
% cdata = imfrmtReformatCellDataSize(cdata, icellfrmt, idatafrmt, ocellfrmt, odatafrmt)
%
% description: 
%     takes a image and reformats it from informat into outformat
%     removes all extra dimensions in img by using first index
%     adds singlentons in missing dimensions
%
% input:
%     cdata        cell data array
%     inDataFrmt   (optional) format of input data (imfrmtFormat(cdata{1}))
%     inCellFrmt   (optional) format of input cell (imfrmtFormat(cdata))
%     outDataFrmt  (optional) format of output data
%     outCellFrmt  (optional) format of output cell
%
% output:
%     cdata        reformatted cell data
%
% note:
%     upper case letters are mathematically in postive axis, lower case mathematically inverted axis

if nargin < 2 || isempty(inCellFrmt)
   inCellFrmt = imfrmtFormatCellSize(cellSize);
end
inCellFrmt = imfrmtFormat(inCellFrmt);

if nargin < 3 || isempty(inDataFrmt)
   inDataFrmt = imfrmtFormat(dataSize);
end
inDataFrmt = imfrmtFormat(inDataFrmt);

if nargin < 4 || isempty(outCellFrmt)
   outCellFrmt = 'SUVW';
end
outCellFrmt = imfrmtFormat(outCellFrmt);

if nargin < 5 || isempty(outDataFrmt)
   outDataFrmt = 'XYZCT';
end
outDataFrmt = imfrmtFormat(outDataFrmt);


if length(unique(lower(inCellFrmt))) ~= length(inCellFrmt)
   error('imfrmtReformatCellSize: labels appear more than once in input cell format: %s', inCellFrmt)
end
if length(unique(lower(inDataFrmt))) ~= length(inDataFrmt)
   error('imfrmtReformatCellSize: labels appear more than once in input data format: %s', inDataFrmt)
end

if length(unique(lower(outCellFrmt))) ~= length(outCellFrmt)
   error('imfrmtReformatCellSize: labels appear more than once in output cell format: %s', outCellFrmt)
end
if length(unique(lower(outDataFrmt))) ~= length(outDataFrmt)
   error('imfrmtReformatCellSize: labels appear more than once in output data format: %s', outDataFrmt)
end

% reformat
ifrmt = [inDataFrmt, inCellFrmt];
ofrmt = [outDataFrmt, outCellFrmt];

asize = [dataSize, cellSize];
asize = imfrmtReformatAnySize(asize, ifrmt, ofrmt);

% split into cells and data
dataSize = asize(1:length(outDataFrmt));
cellSize = asize(length(outDataFrmt)+1:end);

end

