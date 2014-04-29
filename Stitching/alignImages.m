function shifts = alignImages(imgs, param)
%
% shifts = alignImages(imgs, param)
%
% desciption:
%    infers the shifts of images in cell array imgs (rough layout given by cell array structure) 
%    that give rise to a good alignment
% 
% input:
%    imgs    images or image reading commands as cell array, the layout of the cell array should reflect 
%            the neighbouring relations between the images
%    param   (optional) parameter struct with entries:
%            .method       'Optimization', 'SequentialShifts', 'RMS'
%                          + parameters for methods align2ImagesByXXX
%                          ('RMS')
%            .alignment    'Single' = perform single alignment using specified method, 
%                          'All' = perfom alignments in each primary direction of existing neighbours
%                          ('All')        
%
% output: 
%    shifts  absolute shifts of images assuming top left (bottom) image in cell array starts a [0,0(,0)] 
%            in pxel units and pixel coordinates
%
% note:
%    parameter to alignment routines may be specified as a cell to indicate differences in the different 
%    dimensions of cell array imgs
%    for alignment = 'All' SequentialShifts with shift/min/max = 0 is recommended 
%    try to align 3d stacks using a 2d array of 3d images as imgs
%
% See also: align2ImagesByOptimization, align2ImagesBySequentialShift, plotAlignedImages

if nargin < 2
   param = struct();
end

if ~iscell(imgs)
   error('alignImages: expecting cell array as input imgs');
end

adim = ndims1(imgs);
if adim > 3 || adim < 1
   error('alignImages: expect 1d - 3d grid of images, found %gd', adim);
end

algn = getParameter(param, 'alignment', 'All');

meth = getParameter(param, 'method', 'SequentialShifts');

if ischar(meth)  
   meth = repmat({meth}, 1, adim);
end

if ~iscellstr(meth) || length(meth) ~= adim   
   error('alignImages: method not valid');
end

for i = 1:adim
   if ~any(strcmp({'Optimization', 'SequentialShifts', 'RMS'}, meth{i}))
      error('alignImages: method %s not Optimization, SequentialShifts or RMS');
   end
end

% these are the paramter to propagate, check if they are specified as cell arrays
parnames = {'overlap.max', 'overlap.min', 'shift.max', 'shift.min', 'optimizer', 'metric'};
par = repmat({struct()}, 1, adim);

for p = 1:length(parnames)
   parval = getParameter(param, parnames{p}, []);
   if iscell(parval) 
      if length(parval) ~= adim
         erro('alignImages: parameter %s cell array has invalid size',  parnames{i})
      end
      for d = 1:adim
         par{d} = setParameter(par{d}, parnames{p},  parval{d});
      end
   else
      for d = 1:adim
         par{d} = setParameter(par{d}, parnames{p},  parval);
      end
   end
end

%dirs = {'lr', 'tb', 'du'};
dirs = {'tb', 'lr', 'du'};
for d = 1:adim
   par{d} = setParameter(par{1}, 'direction', dirs{d});
end


%note: it might make sense to align form the center image
%      for simplicity we here start at left bottom (down) corner

iread = iscellstr(imgs);
ics = size(imgs);
ics = padright(ics, 3, 1);
shifts = repmat({zeros(1, adim)}, ics);

funs = cell(1,adim);
for a = 1:adim
   funs{a} = str2func(['align2ImagesBy' meth{a}]);
end

imgsc = imgs;

switch algn
   case 'Single' % align in single direction only 
      
      for k = 1:ics(3)
         for j = 1:ics(2)
            for i = 1:ics(1)
               if i == 1
                  if j == 1
                     if k == 1 
                        % i == j == k == 1 
                        if iread
                           % load initial image
                           imgsc{i,j,k} = eval(imgs{i,j,k});
                        end
                     
                     else % i == j == 1, k > 1
                        if iread
                           % clear memory
                           imgsc{ics(1), ics(2), k-1} = imgs{ics(1), ics(2), k-1};
                           
                           %load new image pair
                           imgsc{i,j,k} = eval(imgs{i,j,k});
                           imgsc{i,j,k - 1} = eval(imgs{i,j,k - 1});
                        end
                        
                        shifts{i, j, k} = funs{3}(imgsc{i, j, k-1}, imgsc{i, j, k}, par{3}) + shifts{i, j, k-1};
                        
                        if iread 
                           % clear memory
                           imgsc{i, j, k-1} = imgs{i, j, k-1};
                        end                       
                     end
       
                  else % i == 1, j > 1,  k >= 1;
                     if iread  
                        % clear memory
                        imgsc{ics(1), j-1, k} = imgs{ics(1), j-1, k};
                        
                        %load new image pair
                        imgsc{i,j,k} = eval(imgs{i,j,k});
                        imgsc{i,j-1,k} = eval(imgs{i,j-1,k});
                     end

                     shifts{i, j, k} = funs{2}(imgsc{i, j-1, k}, imgsc{i, j, k}, par{2}) + shifts{i, j-1, k};
                     
                     if iread 
                        % clear memory
                        imgsc{i, j-1, k} = imgs{i, j-1, k};
                     end
                     
                  end
                  
               else % i > 1,  j > 1,  k >= 1; 
                  if iread
                     %load new image
                     imgsc{i,j,k} = eval(imgs{i,j,k});
                  end

                  shifts{i, j, k} = funs{1}(imgsc{i-1, j, k}, imgsc{i, j, k}, par{1}) + shifts{i-1, j, k};
                  
                  if iread   
                     % clear memory
                     imgsc{i-1, j, k} = imgs{i-1, j, k};
                  end
               end
            end
         end
      end
 
      
   case 'All' % algin in all primary directions with corresponding neighbours
      
      for k = 1:ics(3)
         for j = 1:ics(2)
            for i = 1:ics(1)
               if i == 1
                  if j == 1
                     if k == 1 
                        % i == j == k == 1 
                        if iread
                           % load initial image
                           imgsc{i,j,k} = eval(imgs{i,j,k});
                        end
                     
                     else % i == j == 1, k > 1 -> align in k only
                        
                        % align in k
                        if iread
                           % clear memory
                           imgsc{ics(1), ics(2), k-1} = imgs{ics(1), ics(2), k-1};
                           
                           %load new image pair
                           imgsc{i,j,k} = eval(imgs{i,j,k});
                           imgsc{i,j,k - 1} = eval(imgs{i,j,k - 1});
                        end
                        
                        shifts{i, j, k} = funs{3}(imgsc{i, j, k-1}, imgsc{i, j, k}, par{3}) + shifts{i, j, k-1};
                        
                        if iread 
                           % clear memory
                           imgsc{i, j, k-1} = imgs{i, j, k-1};
                        end                       
                     end
       
                  else % i == 1, j > 1,  k >= 1
                     
                     if k == 1 % i == 1, j > 1,  k == 1 -> align in j only
                     
                        % align in j
                        if iread
                           % clear memory
                           imgsc{ics(1), j-1, k} = imgs{ics(1), j-1, k};
                           
                           %load new image pair
                           imgsc{i,j,k} = eval(imgs{i,j,k});
                           imgsc{i,j-1,k} = eval(imgs{i,j-1,k});
                        end
                        
                        shifts{i, j, k} = funs{2}(imgsc{i, j-1, k}, imgsc{i, j, k}, par{2}) + shifts{i, j-1, k};
                        
                        if iread
                           % clear memory
                           imgsc{i, j-1, k} = imgs{i, j-1, k};
                        end
                        
                     else  % i == 1, j > 1,  k > 1 -> align in j and k
                        
                        % align in k
                        if iread
                           % clear memory
                           imgsc{ics(1), j-1, k} = imgs{ics(1), j-1, k};
                           
                           %load new image pair
                           imgsc{i,j,k} = eval(imgs{i,j,k});
                           imgsc{i,j,k - 1} = eval(imgs{i,j,k - 1});
                        end

                        shifts{i, j, k} = funs{3}(imgsc{i, j, k-1}, imgsc{i, j, k}, par{3}) + shifts{i, j, k-1};
                        
                        % align in j
                        if iread 
                           % clear memory
                           imgsc{i, j, k-1} = imgs{i, j, k-1};
                           
                           % load new image
                           imgsc{i,j-1,k} = eval(imgs{i,j-1,k});
                        end
                        
                        sh = funs{2}(imgsc{i, j-1, k}, imgsc{i, j, k}, par{2});
                        shifts{i, j, k}(2) = sh(2) + shifts{i, j-1, k}(2);
                        
                        if iread
                           % clear memory
                           imgsc{i, j-1, k} = imgs{i, j-1, k};
                        end
                        
                     end

                  end
                  
               else % i > 1,  j >= 1,  k >= 1; 
                  
                  if j == 1 % i > 1,  j == 1,  k >= 1;

                     if k == 1 % i > 1,  j == 1,  k == 1 -> align in i

                        % align in i
                        if iread
                           %load new image
                           imgsc{i,j,k} = eval(imgs{i,j,k});
                        end
                        
                        shifts{i, j, k} = funs{1}(imgsc{i-1, j, k}, imgsc{i, j, k}, par{1}) + shifts{i-1, j, k};
                        
                        if iread
                           % clear memory
                           imgsc{i-1, j, k} = imgs{i-1, j, k};
                        end
                     
                     else % i > 1,  j == 1,  k > 1 -> align in i and k
                        
                        % align in i
                        if iread
                           %load new image
                           imgsc{i,j,k} = eval(imgs{i,j,k});
                        end
                        
                        shifts{i, j, k} = funs{1}(imgsc{i-1, j, k}, imgsc{i, j, k}, par{1}) + shifts{i-1, j, k};
                        
                        % align in k
                        if iread
                           % clear memory
                           imgsc{i-1, j, k} = imgs{i-1, j, k};
                           
                           %load new image
                           imgsc{i,j,k - 1} = eval(imgs{i,j,k - 1});
                        end

                        sh = funs{3}(imgsc{i, j, k-1}, imgsc{i, j, k}, par{3});
                        shifts{i, j, k}(3) = sh(3) + shifts{i, j, k-1}(3);
                     end
                     
                  else  % i > 1,  j > 1,  k >= 1;
                     
                     if k == 1 % i > 1,  j > 1,  k == 1 -> align in i and j
                        % align in i
                        if iread
                           %load new image
                           imgsc{i,j,k} = eval(imgs{i,j,k});
                        end
                        
                        shifts{i, j, k} = funs{1}(imgsc{i-1, j, k}, imgsc{i, j, k}, par{1}) + shifts{i-1, j, k};
                        
                        % align in j
                        if iread
                           % clear memory
                           imgsc{i-1, j, k} = imgs{i-1, j, k};
                           
                           %load new image
                           imgsc{i,j-1,k} = eval(imgs{i,j-1,k});
                        end
                        
                        sh = funs{2}(imgsc{i, j-1, k}, imgsc{i, j, k}, par{2});
                        shifts{i, j, k}(2) = sh(2) + shifts{i, j-1, k}(2);
                     
                        if iread
                           % clear memory
                           imgsc{i, j-1, k} = imgs{i, j-1, k};
                        end
  
                     else % i > 1,  j > 1,  k > 1 -> align in i, j and k
                     
                        % align in i
                        if iread
                           %load new image
                           imgsc{i,j,k} = eval(imgs{i,j,k});
                        end
                        
                        shifts{i, j, k} = funs{1}(imgsc{i-1, j, k}, imgsc{i, j, k}, par{1}) + shifts{i-1, j, k};
                        
                        % align in j
                        if iread
                           % clear memory
                           imgsc{i-1, j, k} = imgs{i-1, j, k};
                           
                           %load new image
                           imgsc{i,j-1,k} = eval(imgs{i,j-1,k});
                        end
                        
                        sh = funs{2}(imgsc{i, j-1, k}, imgsc{i, j, k}, par{2});
                        shifts{i, j, k}(2) = sh(2) + shifts{i, j-1, k}(2);
                        
                        
                        % align in k
                        if iread
                           % clear memory
                           imgsc{i, j-1, k} = imgs{i, j-1, k};
                           
                           %load new image
                           imgsc{i,j,k-1} = eval(imgs{i,j,k-1});
                        end
                        
                        sh = funs{3}(imgsc{i, j, k-1}, imgsc{i, j, k}, par{3});
                        shifts{i, j, k}(3) = sh(3) + shifts{i, j, k-1}(3);
                        
                        if iread
                           % clear memory
                           imgsc{i, j, k-1} = imgs{i, j, k-1};
                        end
                     end
                  end
               end
            end
         end
      end
end

end

