function rawSize = imfrmtReshapeInverseSizeOnce(dataSize, dataFrmt, rawFrmt, reshapeFrom, reshapeTo, reshapeSize)
% helper for imfrmtReshape

if length(reshapeFrom) ~= 1
   error('imfrmtReshapeInverseSize: inverse reshaping of ranges for more that a single raw dimension not possible!')
end

if length(reshapeSize) ~=length(reshapeTo)
   error('imfrmtReshapeInverseSize: reshape format %s inconsistent with reshape size %s', reshapeTo, var2char(reshapeSize));
end

posRawFormatFrom = imfrmtPosition(rawFrmt, reshapeFrom);

if length(posRawFormatFrom) ~= length(reshapeFrom)
   error('imfrmtReshapeInverseSize: raw format %s inconsistent with reshape format %s', rawFrmt, reshapeFrom);
end


% output size

rawSize = imfrmtReformatSize(dataSize, dataFrmt, rawFrmt);
rawSize(posRawFormatFrom) = prod(reshapeSize);


end
