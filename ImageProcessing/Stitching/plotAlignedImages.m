function varargout = plotAlignedImages(imgs, shifts, varargin)
%
% img = plotAlignedImages(imgs, shifts)
%
% description:
%    plots aligned images using the shifts and colors conveniently for 
%    inspection of the result
%
% input: 
%    imgs    images
%    shifts  shifts

param = parseParameter(varargin{:});


%size(imgs)
%size(shifts)

if ~iscell(imgs) || ~iscell(shifts) || numel(imgs) ~= numel(shifts)
   error('plotAlignedImages: inconsistent input');
end

% determine image size and absolute shifts w.r.t composed image

imgsizes = cellfun(@size, imgs, 'UniformOutput', false);

[ashifts, asize] = absoluteShiftsAndSize(shifts, imgsizes);

%colors that add up to one
% 1d -> 2 colors, 2d -> 4 colors, 3d -> 8 colors

dim = ndimsd(imgs);
coldim = getParameter(param, 'colors', dim);

switch coldim
   case 1
      cols = {[0, 1, 0], [1, 0, 1]};
   case 2
      cols = {[0.5, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0], [0, 0, 0.5]};
   case 3
      cols =  {...
         [0.25, 0,     0   ],... 
         [0,    0.25,  0   ],...
         [0,    0,     0.25],...
         [0.25  0.25,  0   ],...
         [0,    0.25,  0.25],... 
         [0.25  0,     0.25],... 
         [0.125,0.25,  0   ],...
         [0.125,0,     0.25]};
   otherwise
      cols = {[0.5, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0], [0, 0, 0.5]};
   %   cols = {[0, 1, 0], [1, 0, 1]};
end

si = padright(size(imgs), 3, 1);
cl = length(cols);
img = imgray2color(zeros(asize), 'white');



% switch class(imgs{1})
%    case 'unit32'
%       nrm = 2^32;
%    case 'uint16'
%       nrm = 2^16;
%    case 'uint8'
%       nrm = 2^8;
%    otherwise
%       nrm = 1.0;
% end

imax = cellfun(@(x) max(x(:)), imgs);
nrm = double(max(imax(:)));

cs = 0;
ci = 0;
for k = 1:si(3)
   for j = 1:si(2)
      for i = 1:si(1)
         if i==1
            if j == 1
               cs = mod((k-1) * 4 , 8);
            else
               if si(1) == 1
                  cs = mod(cs + 1, cl);
               else
                  cs = mod(cs + 2, cl);
               end
            end
            ci = 0;
         end
         imga = imgray2color(double(imgs{i, j, k})/nrm, cols{mod(cs + ci, cl)+1});
         imgb = imextract(img, [1 + ashifts{i, j, k}, 1],  [ashifts{i, j, k} + imgsizes{i,j,k}, 3]);
         img  = imreplace(img, imga + imgb, [1 + ashifts{i, j, k},1]);
         ci = ci + 1;
      end
   end
end

implot(img);

if nargout > 0
   varargout{1} = img;
end

end