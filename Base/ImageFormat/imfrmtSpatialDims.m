function d = imfrmtSpatialDims(frmt)
%
% d = imfrmtSpatialDims(frmt)
% d = imfrmtSpatialDims(img)
%
% description: spatial image dimensions from frmt or image

if ~ischar(frmt)
   frmt = imformat(frmt);
end

frmt = lower(frmt);

d = sum(arrayfun(@(x) any(x=='xyz'), frmt));

end