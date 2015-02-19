function si = sized(A, d)
%
% si = sized(A)
%
% description: 
%     calculates size of A adding traling 1s up to dim d
%
% input:
%     A    array
%     d    dim
%
% output:
%     si   size of A with potentialy trailing 1s

si = size(A);
si = padright(si, d, 1);

end
