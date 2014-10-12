function cdata = imfrmtDataToCellData(cdata, dataDims, cellDims) %#ok<INUSD>
%
% cdata = imfrmtDataToCellData(cdata, cellDims)
%
% description: 
%     splits data array to cell array using the last cellDims dimensions
%
% input:
%     cdata      cell data array
%
% output:
%     cdata      reformatted cell data

if ischar(dataDims)
   dataDims = length(dataDims);
end
% if ischar(cellDims)
%    cellDims = length(cellDims);
% end


% split into cells and data
cdata = num2cell(cdata, 1:dataDims);

% remove trailing singleton cell dims
csize = size(cdata);
csize = csize((dataDims+1):end);
if length(csize) <=1
   csize = [csize, 1, 1];
end
cdata = reshape(cdata, csize);

end

