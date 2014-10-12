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

frmt = imfrmtFormat(frmt);
reffrmt = imfrmtFormat(reffrmt);

frmt = lower(frmt);
reffrmt = lower(reffrmt);

[posref, posfrmt] = ismember(reffrmt, frmt);
posfrmt = posfrmt(posref);
posref = find(posref);

end

