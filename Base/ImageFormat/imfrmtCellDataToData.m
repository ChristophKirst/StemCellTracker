function cdata = imfrmtCellDataToData(cdata, dataFrmt, cellFrmt)
%
% cdata = imfrmtCellDataToData(cdata, dataFrmt, cellFrmt)
%
% description: 
%     conversts an c-dim cell array of d-dim data arrays to a c+d dim data array
%     the format will be [datafrmt, cellfrmt]
%
% input:
%     cdata      cell data array
%
% output:
%     cdata      reformatted cell data

% convert data to a big data array

if nargin < 2
   ddim = ndims1(cdata{1});
else
   ddim = length(dataFrmt);
end

if nargin < 3
   cdim = ndims1(cdata);
else
   cdim = length(cellFrmt);
end

csize = size(cdata);
csize = padright(csize, cdim, 1);

cc = cdata;
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

cdata = cc{1};

end

