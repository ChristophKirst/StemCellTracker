function per = impqlformat2permute(frmt1, frmt2)
%
% per = impqlformat2permute(frmt1, frmt2)
%
% description:
%  determines permutation from frmt1 to frmt2 as obtained in imqplpermute, imuvwpermute
%
% input:
%     frmt*   the image/cell formats as chars
%
% output:
%     the pemutation uch that permute(img,per) == impqlreshape(img, frmt1, frmt2)

frmt1 = impqlformat2format(frmt1);
frmt2 = impqlformat2format(frmt2); 

[m,per] = ismember(frmt2, frmt1);
per = per(m);

end