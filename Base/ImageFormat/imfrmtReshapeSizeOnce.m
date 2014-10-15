function outSize = imfrmtReshapeSizeOnce(inSize, inFrmt, outFrmt, reshapeFrom, reshapeTo, reshapeSize)
% helper for imfrmtReshape

if length(reshapeSize) ~=length(reshapeTo)
   error('imfrmtReshape: reshape format %s inconsistent with reshape size %s', reshapeTo, var2char(reshapeSize));
end

posInFormatFrom = imfrmtPosition(inFrmt, reshapeFrom);

if length(posInFormatFrom) ~= length(reshapeFrom)
   error('imfrmtReshape: format %s inconsistent with reshape format %s', inFrmt, reshapeFrom);
end

fromSize = inSize(posInFormatFrom);

if prod(fromSize) ~= prod(reshapeSize)  
   error('imfrmtReshape: size %s inconsistent with reshape size %s', var2char(fromSize), var2char(reshapeSize));  
end

% non reshape format in input format

inFrmtNoReshape = inFrmt;
inFrmtNoReshape(posInFormatFrom) = [];

% reshaped formats in out fromat

[posOutFormatTo,posReshapeTo] = imfrmtPosition(outFrmt, reshapeTo);

if length(posOutFormatTo) ~= length(reshapeTo)
   error('imfrmtReshape: output format %s inconsistent with reshape format %s', outFrmt, reshapeTo);
end

posOutFormatTo = posOutFormatTo(posReshapeTo);

% non reshape format in output format

outFrmtNoReshape = outFrmt;
outFrmtNoReshape(posOutFormatTo) = [];

%if ~strcmp(sort(lower(outFrmtNoReshape)), sort(lower(inFrmtNoReshape)))
%   error('imfrmtReshape: output format %s not consistent with input format %s for reshape from %s to %s',   outFrmtNoReshape, inFrmtNoReshape, reshapeFrom, reshapeTo);
%end

% output size

outSize = imfrmtReformatSize(inSize, inFrmt, outFrmt);
outSize(posOutFormatTo) = reshapeSize;

end
