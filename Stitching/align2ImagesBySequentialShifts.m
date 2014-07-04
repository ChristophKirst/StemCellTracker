function [shift, img] = align2ImagesBySequentialShifts(img1, img2, varargin)
%
% img = align2ImagesBySequentialShifts(img1, img2, param)
%
% description:
%     aligns images img1, img2 optimizing shifts in the directions sequentially
%
% input:
%     img1, img1   images to be alligned, assumed to have same orientation
%     param        (optional)   struct with entries              
%                  .overlap.max maximal overlap of images in primary direction (150)
%                  .overlap.min minimal overlap if images in primary direction (50)
%                  .shift.max   maximal shift in secondary directions ([10, 10])
%                  .shift.min   minimal shift in secondary directions ([-10, -10])
%                  .direction   primary directions: 'lr', 'rl', 'bt', 'tb' 'du', 'ud' 
%                               ({'lr'} = single primary shift assuming img1 left, img2 right)
%                               1, -1, 2, -2, 3, -3 the dimension with sign giving the direction, e.g. 'lr' = 1      
%
% output:
%     shift        the shift between origin of img1 to orgigin img2 in pixel coordinates and pixel units
%     img          (optional) final imaged composed using stitch2ImagesByOverwrite
%
% note: 
%     secondary shifts are 1 in case the primary shift is not 1 or -1 otherwise 2, then the remaining 3rd dimension if image is 3d
%     put shift.min/max to [0,0] to ignore certain secondary directions
%
% See also:  align2ImagesByOptimization, stitch2ImagesByMean, stich2ImagesByWatershed, alignImages, stitchImages, imregister

param = parseParameter(varargin{:});

si1 = size(img1);
si2 = size(img2);

dim = length(si1);
if dim ~= length(si2)
   error('alignImages: images dimension mistmatch!');
end

maxovl = getParameter(param, {'overlap', 'max'}, 150);
if isempty(maxovl) 
   maxovl = 150;
end
minovl = getParameter(param, {'overlap', 'min'}, 50);
if isempty(minovl) 
   minovl = 50;
end


maxsh = getParameter(param, {'shift', 'max'}, ones(1,dim) * 10);
if isempty(maxsh) 
   maxsh = ones(1, dim) * 10;
end
minsh = getParameter(param, {'shift', 'min'}, -ones(1,dim) * 10);
if isempty(minsh) 
   minsh = - ones(1, dim) * 10;
end


direct = getParameter(param, 'direction', 'lr');
if isempty(direct)
   direct = 'lr';
end

if ischar(direct)
   switch dim
      case 2
         if  ~any(strcmp({'lr', 'rl', 'tb', 'bt'}, direct))
            error('align2Images: direction %s not valid');
         end
      case 3
         if  ~any(strcmp({'lr', 'rl', 'tb', 'bt', 'ud', 'du'}, direct))
            error('align2Images: direction %s not valid');
         end
   end
   
   switch direct
      case 'lr'
         direct = 1;
      case 'rl'
         direct = -1;
      case 'bt'
         direct = 2;
      case 'tb'
         direct = -2;
      case 'du'
         direct = 3;
      case 'ud'
         direct = -3;    
   end
   
elseif isnumeric(direct) && length(direct) > 1
   p = find(direct);
   if length(p) > 1
      error('align2Images: direction not valid');
   end
   direct = p * sign(direct(p));
end

if ~isscalar(direct) ||  abs(direct) < 1 || abs(direct) > dim
   error('align2Images: direction not valid'); 
end

shift = zeros(1,dim);

% align in primary direction
minsi = min(si1, si2);
cut1 = cell(1,dim);
cut2 = cut1;
for i = 1:dim
   if i ~= abs(direct)
      cut1{i} = 1:minsi(i);
      cut2{i} = cut1{i};
   else
      cut1{i} = 1:si1(i);
      cut2{i} = 1:si2(i);
   end
end

img1c = img1(cut1{:});
img2c = img2(cut2{:});

% figure(13)
% clf
% imsubplot(1,2,1)
% implot(img1c)
% imsubplot(1,2,2);
% implot(img2c)

bestovl = minovl;
bestdelta = inf;

for ovl = minovl:maxovl
   [img1b, img2b] = extract_overlap(img1c, img2c, ovl, direct, dim);
   
%    size(img1b)
%    size(img2b)

   delta = sum(sum(abs(img1b-img2b))) / ovl / mean(img1b(:)) / mean(img2b(:));
    if delta < bestdelta
        bestovl = ovl;
        bestdelta = delta;
    end
end

if direct > 0
   shift(direct) = size(img1, direct) - bestovl;
else
   shift(-direct) = -size(img2, -direct) + bestovl;
end

% align in secondary direction

[img1c, img2c] = extract_overlap(img1, img2, bestovl, direct, dim);

if abs(direct) == 1
   d = 2;
else
   d = 1;
end

if dim == 3
   dr = setdiff([1 2 3], [direct, d]);
   minsi = min(si1(dr), si2(dr));
    
   cut1 = cell(1,dim);
   cut2 = cut1;
   for i = 1:dim
      if i ~= abs(direct)
         cut1{i} = 1:bestovl;
         cut2{i} = cut1{i};
      elseif i == dr
         cut1{i} = 1:minsi;
         cut2{i} = cut1{i};
      else
         cut1{i} = 1:si1(i);
         cut2{i} = 1:si2(i);
      end
   end
   
   img1d = img1c(cut1{:});
   img2d = img2c(cut2{:});
else
   img1d = img1c;
   img2d = img2c;
end

% figure(14)
% clf
% imsubplot(1,2,1)
% implot(img1d)
% imsubplot(1,2,2);
% implot(img2d)
% 
% size(img1d)
% size(img2d)


bestsh = minsh(1);
bestovl = inf;
bestdelta = inf;

sid1 = size(img1d, d);
sid2 = size(img2d, d);

for sh = minsh(1):maxsh(1)
   
   if sh > 0
      ovl = min(sid1-sh, sid2);  
   else
      ovl = min(sid2+sh, sid1);  
   end
   
   if ovl >= 0
      [img1b, img2b] = extract_overlap(img1d, img2d, ovl, signp(sh) * d, dim);
      
      delta = sum(sum(abs(img1b-img2b))) / ovl / mean2(img1b) / mean2(img2b);
      if delta < bestdelta
         bestsh = sh;
         bestovl = ovl;
         bestdelta = delta;
      end
   end
end

shift(d) = bestsh;

if dim > 2

   [img1c, img2c] = extract_overlap(img1c, img2c, bestovl, signp(bestsh) * d, dim);
   
%    figure(15)
%    clf
%    imsubplot(1,2,1)
%    implot(img1c)
%    imsubplot(1,2,2);
%    implot(img2c)

   bestsh = minsh(2);
   bestdelta = inf;

   for sh = minsh(2):maxsh(2)
      
      if sh > 0
         ovl = si(dr) - sh;
      else
         ovl = si(dr) - sh;
      end
      
      if ovl >= 0
         [img1b, img2b] = extract_overlap(img1c, img2c, ovl, signp(sh) * dr, dim);
         
         delta = sum(sum(abs(img1b-img2b))) / ovl / mean2(img1b) / mean2(img2b);
         if delta < bestdelta
            bestsh = sh;
            bestdelta = delta;
         end
      end
   end
   
   shift(dr) = bestsh;
end
   
   
if nargout > 1
   img = stitch2ImagesByOverwrite(img1, img2, shift);
end

end


function [img1b, img2b] = extract_overlap(img1, img2, ovl, direct, dim)
% extract the borders
switch dim
   case 2
      switch direct
         case 1
            img1b = img1(end-ovl+1:end, :);
            img2b = img2(1:ovl, :);
         case -1
            img2b = img2(end-ovl+1:end,:);
            img1b = img1(1:ovl,:);
         case 2
            img1b = img1(:, end-ovl+1:end);
            img2b = img2(:, 1:ovl);
         case -2
            img2b = img2(:, end-ovl+1:end);
            img1b = img1(:, 1:ovl);
         otherwise
            error('align2Images: direction not valid!');
      end
      
   case 3
      switch direct
         case 1
            img1b = img1(end-ovl+1:end, :, :);
            img2b = img2(1:ovl, :, :);
         case -1
            img2b = img2(end-ovl+1:end,:, :);
            img1b = img1(1:ovl,:, :);
         case 2
            img1b = img1(:, end-ovl+1:end, :);
            img2b = img2(:, 1:ovl, :);
         case -2
            img2b = img2(:, end-ovl+1:end, :);
            img1b = img1(:, 1:ovl+1, :);
         case 3
            img1b = img1(:, :, end-ovl+1:end, :);
            img2b = img2(:, 1:ovl+1, :);
         case -3
            img2b = img2(:, end-ovl+1:end, :);
            img1b = img1(:, 1:ovl+1, :);
         otherwise
            error('align2Images: direction not valid!');
      end
end   
end

         


