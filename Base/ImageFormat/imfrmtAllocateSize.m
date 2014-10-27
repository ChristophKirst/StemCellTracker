function si = imfrmtAllocateSize(si)
%  
% si = imfrmtAllocateSize(si)
%
% make si to accomply with the matlab allocation of cells / zeros
%

if length(si) < 2
    si = [si, 1];
end

end