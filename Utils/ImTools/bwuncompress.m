function [mask, isize] = bwuncompress(mm, type, format)
%
% mm = bwuncompress(mm, type, format)
%
% description:
%   undo compression of bwimage
%
% input:
%   mm      compressed image as obtained from bwcompress
%   type    (optional) 'full','edge' ([] = 'full')
%   format  (optional) 'mask', 'list'='PixelIdxList' ('mask')
%
% output:
%   mask    the bw mask or pixel list
%   isize   (optional) image size
%
% note:
%     The 'edge' option does not completely fill in the vertical edges, need
%     take full mask, erode and subtract.
%
% See also: bwcompress, bwpack

if nargin < 2
   type = 'full';
end
if isempty(type)
   type = 'full';
end
type = lower(type);
if ~any(strcmp({'full', 'edge'}, type))
   error('bwuncompress: type is not full or edge but %s', type);
end

if nargin < 3
   format = 'mask';
end
if isempty(format)
   format = 'mask';
end
format = lower(format);
if strcmp(format, 'pixelidxlist') 
   format = 'list';
end
if ~any(strcmp({'list', 'mask'}, format))
   error('bwuncompress: format is not mask, list or PixelIdxList but %s', format);
end

dim = mm(1);
isize = mm(2:dim+1);

switch format
   case 'mask'
      
      mask = false(isize);
      switch type
         case 'full'
            for i = dim+2:2:length(mm)
               mask(mm(i):mm(i+1)) = 1;
            end
         case 'edge'
            if strcmp(type, 'edge')
               mask(mm(dim+2:end)) = 1;
            end
      end
 
   case 'list'
      switch type
         case 'full'
            dl = mm(dim+3:2:end) - mm(dim+2:2:end) + 1;  %number of positions in each row
            nl = sum(dl);
            mask = zeros(1,nl);
            dl = [0 cumsum(dl)];
            ioff = dim + 1;
            for i = 1:length(dl)-1         
               mask(dl(i)+1:dl(i+1)) = mm(2*i-1 + ioff):mm(2*i + ioff);
            end
         case 'edge'
            mask = mm(dim+2:end);
      end  
end

end