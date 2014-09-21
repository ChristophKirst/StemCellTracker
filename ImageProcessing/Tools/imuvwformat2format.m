function frmt = imuvwformat2format(frmt)
%
% frmt = imuvwformat2format(frmt)
%
% description:
%     changes the image format by replacing x,y,z -> u,v,w
% 

frmt(frmt == 'x') = 'u';
frmt(frmt == 'y') = 'v';
frmt(frmt == 'z') = 'w';

end