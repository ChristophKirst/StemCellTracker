function per = imfrmtPermuation(infrmt, outfrmt)
%
% per = imfrmtPermuation(infrmt, outfrmt)
%
% description:
%  determines permutation from infrmt to outfrmt (up to axis inversion)
%
% input:
%     *frmt   the image/cell formats as chars
%
% output:
%     the pemutation such that permute(img,per) == imfrmtReformat(img, frmt1, frmt2) up to axis inversions

if ~iscellstr(infrmt)
   infrmt = imfrmtFormat(infrmt);
end
if ~iscellstr(outfrmt)
   outfrmt = imfrmtFormat(outfrmt); 
end

[m,per] = ismember(lower(outfrmt), lower(infrmt));
per = per(m);

end