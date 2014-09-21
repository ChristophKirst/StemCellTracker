function frmt = impqlformat2format(frmt)
%
% frmt = impqlformat2format(frmt)
%
% description:
%     changes the image format by replacing x,y,z -> p,q,l
% 

frmt(frmt == 'x') = 'p';
frmt(frmt == 'y') = 'q';
frmt(frmt == 'z') = 'l';

end