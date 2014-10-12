function cellData = imfrmtReformatCellData(cellData, inDataFrmt, inCellFrmt, outDataFrmt, outCellFrmt)
%
% cellData = imfrmtReformatCellData(cellData, inDataFrmt, inCellFrmt, outDataFrmt, outCellFrmt)
%
% description: 
%     takes a image and reformats it from informat into outformat
%     removes all extra dimensions in img by using first index
%     adds singlentons in missing dimensions
%
% input:
%     cellData     cell data array
%     inDataFrmt   (optional) format of input data (imfrmtFormat(cdata{1}))
%     inCellFrmt   (optional) format of input cell (imfrmtFormat(cdata))
%     outDataFrmt  (optional) format of output data
%     outCellFrmt  (optional) format of output cell
%
% output:
%     cellData   reformatted cell data
%
% note:
%     upper case letters are mathematically in postive axis, lower case mathematically inverted axis

if nargin < 2 || isempty(inCellFrmt)
   inCellFrmt = imfrmtFormat(cellData);
end
inCellFrmt = imfrmtFormat(inCellFrmt);

if nargin < 3 || isempty(inDataFrmt)
   inDataFrmt = imfrmtFormat(cellData{1});
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
   error('imfrmtReformatCell: labels appear more than once in input cell format: %s', inCellFrmt)
end
if length(unique(lower(inDataFrmt))) ~= length(inDataFrmt)
   error('imfrmtReformatCell: labels appear more than once in input data format: %s', inDataFrmt)
end

if length(unique(lower(outCellFrmt))) ~= length(outCellFrmt)
   error('imfrmtReformatCell: labels appear more than once in output cell format: %s', outCellFrmt)
end
if length(unique(lower(outDataFrmt))) ~= length(outDataFrmt)
   error('imfrmtReformatCell: labels appear more than once in output data format: %s', outDataFrmt)
end


% convert data to a big data array
cdim = length(inCellFrmt);
ddim = length(inDataFrmt);
csize = imfrmtSize(cellData, inCellFrmt);

cc = cellData;
for cd = 1:cdim
   
   cs = csize(cd+1:end);
   nc = prod(cs);
   
   asgn = num2cell(ones(1, cdim-cd+1));
   asgn{1} = ':';
   
   if length(cs) == 1
      ccn = cell(cs,1);
   else
      ccn = cell(cs);
   end

   for i = 1:nc
      sub = imind2sub(cs, i);
      asgn(2:end) = num2cell(sub);
      cca = cc(asgn{:});
      ccn{i} = cat(ddim + cd, cca{:});
   end

   cc = ccn;
end

cc = cc{1};

% reformat
ifrmt = [inDataFrmt, inCellFrmt];
ofrmt = [outDataFrmt, outCellFrmt];

cc = imfrmtReformatAny(cc, ifrmt, ofrmt);

% split into cells and data

cdim = length(outDataFrmt);
cellData = num2cell(cc, 1:cdim);

% remove trailing singleton cell dims
csize = size(cellData);
csize = csize((cdim+1):end);
if length(csize) <= 1
   csize = [csize, 1,1];
end
cellData = reshape(cellData, csize);

end

