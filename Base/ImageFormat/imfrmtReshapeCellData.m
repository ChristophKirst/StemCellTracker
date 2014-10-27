function outData = imfrmtReshapeCellData(inCellData, inDataFrmt, inCellFrmt, outDataFrmt, outCellFrmt, reshapeFrom, reshapeTo, reshapeSize)
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
%    reshaping from more than two coordinates is possible in principle but
%    usally not required for image data
%    reshapeFrom is thus restricted to be a single character for simplicity


% determine the strucure of the reshaping, pure cell vs. data / cell mixtures
[inPureCellType, outPureCellType, ...
 cellReshapeFrom, cellReshapeTo, cellReshapeSize, ...
 reducedReshapeFrom, reducedReshapeTo, reducedReshapeSize] = ...
    formatType(inDataFrmt, inCellFrmt, outDataFrmt, outCellFrmt, ...
               reshapeFrom, reshapeTo, reshapeSize);
%             
%   inPureCellType    
%   outPureCellType

% move pure cell types to end of cell frmt

inCellFrmtReshape = inCellFrmt(~inPureCellType);
inCellFrmtPure = inCellFrmt(inPureCellType);
inCellFrmtOrdered = [inCellFrmtReshape, inCellFrmtPure];

inCellData = imfrmtReformat(inCellData, inCellFrmt, inCellFrmtOrdered);

inDataSize = imfrmtSize(inCellData{1}, inDataFrmt);
inCellSize = imfrmtSize(inCellData, inCellFrmtOrdered);


% cell format after reshaping

outCellFrmtReshape = outCellFrmt(~outPureCellType);
%outCellFrmtPure = outCellFrmt(outPureCellType);
outCellFrmtTemp = [outCellFrmtReshape, inCellFrmtPure];

% reshape inner mixture data 

n = prod(inCellSize(inPureCellType));

[~, outCellSizeOrdered] = imfrmtReshapeCellDataSize(inDataSize, inCellSize,...
                              inDataFrmt, inCellFrmtOrdered, outDataFrmt, outCellFrmtTemp, reducedReshapeFrom, reducedReshapeTo, reducedReshapeSize);

outCellSizeOrdered = imfrmtAllocateSize(outCellSizeOrdered);
                           
outData = cell(outCellSizeOrdered);
%size(outData)


argsIn = repmat({':'}, 1, sum(~inPureCellType) + 1);
argsOut = repmat({':'}, 1, sum(~outPureCellType) + 1);

fullInFrmt = [inDataFrmt, inCellFrmtReshape];
fullOutFrmt = [outDataFrmt, outCellFrmtReshape];

for i = 1:n
   argsIn{end} = i;
   argsOut{end} = i;

   % convert to full data array
   data = imfrmtCellDataToData(inCellData(argsIn{:}), inDataFrmt, inCellFrmtReshape);

   % reshape
   data =imfrmtReshape(data, fullInFrmt, fullOutFrmt, reducedReshapeFrom, reducedReshapeTo, reducedReshapeSize);

   % convert to cell array
   data = imfrmtDataToCellData(data, outDataFrmt, outCellFrmtReshape);
   
   % assing to Data
   outData(argsOut{:}) = data;
end


% size(outData)
% outCellFrmtTemp
% outCellFrmt
% 
% cellReshapeFrom
% cellReshapeTo
% cellReshapeSize

% reshape cell data

outData = imfrmtReshape(outData, outCellFrmtTemp, outCellFrmt, cellReshapeFrom, cellReshapeTo, cellReshapeSize);

end






function [inCellType, outCellType, cellReshapeFrom, cellReshapeTo, cellReshapeSize, reducedReshapeFrom, reducedReshapeTo, reducedReshapeSize] = formatType(inDataFrmt, inCellFrmt, outDataFrmt, outCellFrmt, reshapeFrom, reshapeTo, reshapeSize)

% classifies the frmt into dims that are data/cell mxiture types and pure cell types 
% for pure cell types simple reshaping is faster
% assumes full data / cell output formats !, i.e. singleton dimensions not removed !

reshapeFrom = lower(reshapeFrom);
reshapeTo   = lower(reshapeTo);

inDataFrmt  = lower(inDataFrmt);
inCellFrmt  = lower(inCellFrmt);

outDataFrmt  = lower(outDataFrmt);
outCellFrmt  = lower(outCellFrmt);

n = length(reshapeFrom);

% set input types

frmtIn = struct;
for i = 1:length(inDataFrmt)
   frmtIn.(inDataFrmt(i)) = 'd';
end
for i = 1:length(inCellFrmt)
   frmtIn.(inCellFrmt(i)) = 'c';
end


% propagate to output frmt
frmtOut = frmtIn;
for r = 1:n
   rf = reshapeFrom{r};
   f = frmtOut.(rf(1));
%    for i = 2:length(rf)  % we assume reshapeFrom to bi single dim !
%       if f ~= frmtOut.(rf(i))
%          f = 'b';
%       end
%    end
%   frmtOut = rmfield(frmtOut, num2cell(rf));

   frmtOut = rmfield(frmtOut, num2cell(rf));
   
   rt = reshapeTo{r};
   for i = 1:length(rt)
      frmtOut.(rt(i)) = f;
   end
end

% chect for type of progagated types
for i = 1:length(outDataFrmt)
   if frmtOut.(outDataFrmt(i)) ~= 'd'
      frmtOut.(outDataFrmt(i)) = 'b';
   end
end
for i = 1:length(outCellFrmt)
   if frmtOut.(outCellFrmt(i)) ~= 'c'
      frmtOut.(outCellFrmt(i)) = 'b';
   end
end


% back propagate types 
for r = n:-1:1
   rt = reshapeTo{r};
   rf = reshapeFrom{r};

   f = frmtOut.(rt(1));
   for i = 2:length(rt)  
      if f ~= frmtOut.(rt(i))
         f = 'b';
      end
   end
   frmtOut = rmfield(frmtOut, num2cell(rt)); 
   frmtOut.(rf) = f;
end

frmtIn = frmtOut;

% final forward propagate for consistent types and reshape froms

cellReshapeFrom = {};
reducedReshapeFrom = {};

cellReshapeTo = {};
reducedReshapeTo = {};

cellReshapeSize = {};
reducedReshapeSize = {};

for r = 1:n
   rt = reshapeTo{r};
   rf = reshapeFrom{r};
   f = frmtOut.(rf);
   frmtOut = rmfield(frmtOut, rf);  
   for i = 1:length(rt)
      frmtOut.(rt(i)) = f;
   end  
   
   if f == 'c'
      cellReshapeFrom{end+1} = reshapeFrom{r}; %#ok<AGROW>
      cellReshapeTo{end+1}   = reshapeTo{r}; %#ok<AGROW>
      cellReshapeSize{end+1} = reshapeSize{r}; %#ok<AGROW>
   else
      reducedReshapeFrom{end+1} = reshapeFrom{r};%#ok<AGROW>
      reducedReshapeTo{end+1}   = reshapeTo{r};%#ok<AGROW>
      reducedReshapeSize{end+1} = reshapeSize{r};%#ok<AGROW>
   end
end


n = length(inCellFrmt);
inCellType = false(1, n, 'logical');
for i = 1:n
   inCellType(i) = frmtIn.(inCellFrmt(i)) == 'c';
end

n = length(outCellFrmt);
outCellType = false(1, n, 'logical');
for i = 1:n
   outCellType(i) = frmtOut.(outCellFrmt(i)) == 'c';
end

   
end



