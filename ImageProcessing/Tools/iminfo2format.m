function format = iminfo2format(iinfo)
%
% format = iminfo2format(iinfo)
%
% description:
%     determines the format of an image given the iinfo struct that is the reuslt of imread_bf_info
% 
% input:
%     iinfo    image info as returned by imread_bf_info
%
% output:
%     format:
%     'pq'   = p x q matrix = 'pq' (2D grayscale)
%     'pqwc'  = p x q x c matrix (2D multi channel )
%     'pql'  = p x q x l matrix (3D)
%     'pqcl' = p x q x c x l matrix (3D multi channel image, matlab ordering)
%     'pqlc' = p x q x l x c matrix (3D multi channel image)
%     'pqlt'  = p x q x l x t matrix (4D grayscale)
%     'pqclt' = p x q x c x l x t matrix (4D multi channel image, matlab ordering)
%     'pqlct' = p x q x l x c x t matrix (4D multi channel image, time last)
%     'pqlct' = p x q x l x t x c matrix (4D multi channel image, channel last)
%     ''  = not supported format
%     p = x pixel coordinate, q = z pixel coordinate, l = z pixel coordinate, c = color, t = time
%
% note:
%     'pql' = 'pqt', 'pqcl' = 'pqct', 'pqlc' = 'pqtc'
%
% See also: imformat






end
      
   




