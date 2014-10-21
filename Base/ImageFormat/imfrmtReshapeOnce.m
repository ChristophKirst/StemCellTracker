function data = imfrmtReshapeOnce(data, inFrmt, outFrmt, reshapeFrom, reshapeTo, reshapeSize)
% helper for imfrmtReshape

% check reshape specs
if length(reshapeSize) ~=length(reshapeTo)
   error('imfrmtReshape: reshape format %s inconsistent with reshape size %s', reshapeTo, var2char(reshapeSize));
end

posInFormatFrom = imfrmtPosition(inFrmt, reshapeFrom);
if length(posInFormatFrom) ~= length(reshapeFrom)
   error('imfrmtReshape: format %s inconsistent with reshape format %s', inFrmt, reshapeFrom);
end

%data
%inFrmt
inSize = imfrmtSize(data, inFrmt);
fromSize = inSize(posInFormatFrom);
if prod(fromSize) ~= prod(reshapeSize)  
   error('imfrmtReshape: size %s inconsistent with reshape size %s', var2char(fromSize), var2char(reshapeSize));  
end

% reshaped formats in out fromat
posOutFormatTo = imfrmtPosition(outFrmt, reshapeTo);
if length(posOutFormatTo) ~= length(reshapeTo)
   error('imfrmtReshape: output format %s inconsistent with reshape format %s', outFrmt, reshapeTo);
end


%  reformat all from formats to single block
reshapeBase = imfrmtRemoveFormat(inFrmt, reshapeFrom);
reshapeInFrmt = [reshapeBase, reshapeFrom];
reshapeOutFrmt =  [reshapeBase, reshapeTo];
data = imfrmtReformat(data, inFrmt, reshapeInFrmt);

% reshape size
rSize = imfrmtSize(data, reshapeBase);
rSize = [rSize, reshapeSize];

data = reshape(data, rSize);

% figure(300);
% implottiling(imfrmtReformatCellData({data}, 'XYS', 'Q', 'XY', 'S'));

data = imfrmtReformat(data, reshapeOutFrmt, outFrmt);

% figure(301);
% implottiling(imfrmtReformatCellData({data}, 'XYS', 'Q', 'XY', 'S'))
% reshapeOutFrmt
% outFrmt


end
