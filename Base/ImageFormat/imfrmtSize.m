function si = imfrmtSize(d, frmt)
%
% si = imfrmtSize(d, frmt)
%
% description:
%   returns size assuming a format frmt

if ischar(frmt)
   n = length(frmt);
else
   n = frmt;
end

si = size(d);
si = padright(si, n, 1);


end