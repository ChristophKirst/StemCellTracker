function d = imsdims(iformat)
%
% d = imsdims(iformat)
%
% description: spatial image dimensions from iformat or image

if ~ischar(iformat)
   iformat = imformat(iformat);
end

d = sum(arrayfun(@(x) any(x=='xyzpql'), iformat));

end