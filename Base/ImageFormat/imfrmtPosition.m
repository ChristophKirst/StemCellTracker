function [posref, posfrmt]  = imfrmtPosition(reffrmt, frmt)
%
% pos = imfrmtPosition(reffrmt, frmt)
%
% description:
%     detremines positions of frmt dimensions in reffrmt, ignores axis inversions
% 
% input:
%     reffrmt  reference format 
%     frmt     format 
% 
% output:
%     pos      positions of the chracters in frmt in reffrmt, 0 if not found 
%
% See also: imfrmtFormat

if ~iscellstr(frmt)
   frmt = imfrmtFormat(frmt);
end
if ~iscellstr(reffrmt)
   reffrmt = imfrmtFormat(reffrmt);
end

frmt = lower(frmt);
reffrmt = lower(reffrmt);

[posref, posfrmt] = ismember(reffrmt, frmt);
posfrmt = posfrmt(posref);
posref = find(posref);

end

