function [ovl1, ovl2] = overlap2ImagesOnGrid(imgs, varargin)
%
% [ovl1, ovl2] = overlap2ImagesOnGrid(imgs, param)
%
% description:
%    returns the sub-images of img1 and img2 that contribute to the potential overlap
%    given parameter overlap.max and alignemnt on gird
%
% intput:
%    imgs        images a pre aligned cell array
%    param       parameter struct with entries
%                .overlap.max
%
% output:
%    ovl1, ovl2   the potential overlapping region extracted form imgs


if numel(imgs) ~= 2
   error('overlap2ImagesOnGrid: expect cell array with two images!');
end

param = parseParameter(varargin{:});

mo = getParameter(param, 'overlap.max', []);

img1 = imgs{1}; img2 = imgs{2};

% orientation

si = size(imgs);
dim = ndims(imgs);
pos = find(si == 2, 1);
idx = repmat({':'}, 1, dim);

% img1 overlap region
si1 = size(img1);
sip  = si1(pos);

if isempty(mo)
   moi = sip;
else
   moi = mo;
end

idx{pos} = sip-moi+1:sip;

ovl1 = img1(idx{:});

% img2 overlap region

si2 = size(img2);
sip  = si2(pos);

if isempty(mo)
   moi = sip;
else
   moi = mo;
end

idx{pos} = 1:moi;

ovl2 = img2(idx{:});

end



