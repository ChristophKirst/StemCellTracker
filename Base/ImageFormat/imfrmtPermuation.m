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

infrmt = imfrmtFormat(infrmt);
outfrmt = imfrmtFormat(outfrmt); 

[m,per] = ismember(outfrmt, infrmt);
per = per(m);

end