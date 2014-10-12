function [extrafrmt, extrapos]  = imfrmtExtraFormats(reffrmt, frmt)
%
% pos = imfrmtExtraFormats(reffrmt, frmt)
%
% description:
%     determines additional formats in frmt w.r.t. the reffrmt, ignores axis inversions
% 
% input:
%     reffrmt   reference format 
%     frmt      format 
% 
% output:
%     extrafrmt extra dimensions
%     extrapos  positions in frmt of extra dimensions
%
% See also: imfrmtFormat

frmt = imfrmtFormat(frmt);
reffrmt = imfrmtFormat(reffrmt);

frmtl = lower(frmt);
reffrmt = lower(reffrmt);

posref = ismember(frmtl, reffrmt);
extrapos = find(posref==0);
extrafrmt = frmt(extrapos);

end

