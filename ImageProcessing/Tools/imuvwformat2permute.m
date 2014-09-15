function per = imuvwformat2permute(frmt1, frmt2)
%
% per = imuvwformat2permute(frmt1, frmt2)
%
% description:
%  determines permutation from frmt1 to frmt2 as obtained in imqplpermute, imuvwpermute
%
% input:
%     frmt*   the image/cell formats as chars
%
% output:
%     the permutation such that permute(img,per) == imuvwreshape(img, frmt1, frmt2)

frmt1(frmt1 == 'x') = 'u';
frmt1(frmt1 == 'y') = 'v';
frmt1(frmt1 == 'z') = 'w';

frmt2(frmt2 == 'x') = 'u';
frmt2(frmt2 == 'y') = 'v';
frmt2(frmt2 == 'z') = 'w';

[m,per] = ismember(frmt2, frmt1);
per = per(m);

end