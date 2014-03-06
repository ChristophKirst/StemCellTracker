function formats = imformatlist()
%
% formats = imformatlist()
%
% description:
%     returns possible image formats
% 
% output:
%     format:
%     'hw'   = h x w matrix = 'hw' (2D grayscale)
%     'hwc'  = h x w x c matrix (2D multi channel )
%     'hwl'  = h x w x l matrix (3D)
%     'hwcl' = h x w x c x l matrix (3D multi channel image, matlab ordering)
%     'hwlc' = h x w x l x c matrix (3D multi channel image)
%     'hwlt'  = h x w x l x t matrix (4D grayscale)
%     'hwclt' = h x w x c x l x t matrix (4D multi channel image, matlab ordering)
%     'hwlct' = h x w x l x c x t matrix (4D multi channel image, time last)
%     'hwlct' = h x w x l x t x c matrix (4D multi channel image, channel last)
%     ''  = not supported format
%     h = height, w = width, l = length, c = color, t = time
%
% note:
%     'hwl' = 'hwt', 'hwcl' = 'hwct', 'hwlc' = 'hwtc'

formats = {'hw', 'hwc', 'hwl', 'hwcl', 'hwlc', 'hwlt', 'hwclt', 'hwlct', 'hwltc'};

end
   




