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

frmt1(frmt1 == 'x') = 'p';
frmt1(frmt1 == 'y') = 'q';
frmt1(frmt1 == 'z') = 'l';

frmt2(frmt2 == 'x') = 'p';
frmt2(frmt2 == 'y') = 'q';
frmt2(frmt2 == 'z') = 'l';

[m,per] = ismember(frmt2, frmt1);
per = per(m);

end