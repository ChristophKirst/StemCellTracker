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

if ~iscell(frmt)
   frmt = imfrmtFormat(frmt);
end
if ~iscell(reffrmt)
   reffrmt = imfrmtFormat(reffrmt);
end

frmtl = lower(frmt);
reffrmt = lower(reffrmt);

posref = ismember(frmtl, reffrmt);
extrapos = find(posref==0);
extrafrmt = frmt(extrapos);

end

