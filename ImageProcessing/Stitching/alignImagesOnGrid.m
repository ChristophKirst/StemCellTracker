function shifts = alignImagesOnGrid(imgs, varargin)
%
% shifts = alignImagesOnGrid(imgs, varargin)
%
% desciption:
%    infers shifts of images in cell array imgs to make them align
%    the fucntion is a wrapper that distributes the tasks to the indibvidual align and optimization methods
% 
% input:
%    imgs    images or image reading commands as cell array, the layout of the cell array usually should reflect 
%            the neighbouring relations between the images
%    param   (optional) parameter struct with entries:
%            .alignment       'Optimization', 'RMS', 'Correlation', 'Hugin' ('RMS')
%                             + parameters for methods align2ImagesOnGridByXXX
%            .method         'sequential' = perform single alignment using specified method, 
%                             'primary'= perform alignments in each primary direction of existing neighbours
%                             'global' = perform alignment globally given the individual pairwise alignments
%                             ('global')
%            .anchor          starting image used if method is 'sequential' or 'primary', shifts are recursivle 
%                             propagted from this image to nearest neighbours
%
% output: 
%    shifts  absolute shifts of images assuming [lower,left(,bottom)] of image in cell array starts a [0,0(,0)] 
%            in pixel units and pixel coordinates
%
% See also: alignImages, plotAlignedImages



alignnames = {'Optimization', 'Correlation', 'RMS', 'Hugin'};
methodnames = {'sequential', 'primary', 'global'};

param = parseParameter(varargin{:});

if ~iscell(imgs)
   error('alignImagesOnGrid: expecting cell array as input imgs');
end

adim = ndims1(imgs);
if adim > 3 || adim < 1
   error('alignImagesOnGrid: expect 1d - 3d grid of images, found %gd', adim);
end

algn = getParameter(param, 'alignment', 'RMS');
meth = getParameter(param, 'method', 'global');

if ~any(strcmp(methodnames, meth))
      error('alignImagesOnGrid: method %s not in %s', meth, var2char(methodnames));
end

if ~any(strcmp(alignnames, algn))
      error('alignImagesOnGrid: alignment %s not in %s', align, var2char(alignnames));
end


fun = str2func(['align2ImagesOnGridBy' algn]);

shifts = {};
switch meth
   case 'sequential'
      shifts = alignBySequentialShifts(imgs, fun, param);
      
   case 'primary'
      shifts = alignByPrimary(imgs, fun, param);
      
   case 'global'
      shifts = alignByGlobal(imgs, fun, param);     
end    


end


% helper

%%%%%%%%%%%%%%%%%%%%%%%%
%%% sequential shifts 

function shifts = alignBySequentialShifts(imgs, fun, param)

si = size(imgs); dim = ndims(imgs);
start = getParameter(param, 'anchor', ceil(si/2));
cstart = num2cell(start);

shifts = cell(si);
shifts{cstart{:}} = zeros(1,dim);
shifts = alignNeighbours(imgs, shifts, start, dim, si, fun, param);

shifts = zero2first(shifts);

end


function shifts = alignNeighbours(imgs, shifts, start, dim, si, fun, param)

   % nearest neighbours
   cstart = num2cell(start); 
   rec = {};
   for d = 1:dim
      if start(d) < si(d)
         neigb = start;
         neigb(d) = start(d) + 1;
         neigb = num2cell(neigb);
         
         if isempty(shifts{neigb{:}})
            pos = ones(1,dim); pos = num2cell(pos);
            clear pair;
            pair{pos{:}} = imgs{cstart{:}}; %#ok<AGROW>
            pos{d} = pos{d} + 1;
            pair{pos{:}} = imgs{neigb{:}}; %#ok<AGROW> 
            
            shifts{neigb{:}} = fun(pair, param) + shifts{cstart{:}};
            rec = {rec{:}, [neigb{:}]}; %#ok<CCAT>
         end
      end
      
      if start(d) > 1
         neigb = start;
         neigb(d) = start(d) - 1;
         neigb = num2cell(neigb);
         
         if isempty(shifts{neigb{:}})
            pos = ones(1,dim); pos = num2cell(pos);
            clear pair;
            pair{pos{:}} = imgs{neigb{:}}; %#ok<AGROW>
            pos{d} = pos{d} + 1;
            pair{pos{:}} = imgs{cstart{:}}; %#ok<AGROW>
   
            shifts{neigb{:}} = - fun(pair, param) + shifts{cstart{:}};
            rec = {rec{:}, [neigb{:}]}; %#ok<CCAT>
         end
      end
   end
   
   % recursion 
   for r = 1:length(rec)
      shifts = alignNeighbours(imgs, shifts, rec{r}, dim, si, fun, param);
   end
   
end


% zero first to first shift
function shifts = zero2first(shifts)

   sh = shifts{1};
   for i = 1:numel(shifts)
      shifts{i} = shifts{i} - sh;
   end

end


%%%%%%%%%%%%%%%%%%%%%%%%
%%% primary direction

function sh = alignByPrimary(imgs, fun, param)

   % determine pairs 
   si = size(imgs); dim = ndims(imgs);
   si = padright(si, 3, 1);

   np = si(1) * si(2) * si(3);
   shifts = cell([np, np, np]);
   
   for x = 1:si(1)
      for y = 1:si(2)
         for z = 1:si(3)
            if x < si(1)
               from = imsub2ind(si, [x,y,z]);
               to = imsub2ind(si, [x+1,y,z]);
               shifts{from,to} = {imgs{x,y,z}; imgs{x+1, y, z}};         
            end
            if y < si(2)
               from = imsub2ind(si, [x,y,z]);
               to = imsub2ind(si, [x,y+1,z]);
               shifts{from, to} = {imgs{x,y,z}, imgs{x, y+1, z}};
            end
            if z < si(3)
               from = imsub2ind(si, [x,y,z]);
               to = imsub2ind(si, [x,y,z+1]);
               c{1,1,1} = imgs{x,y,z}; c{1,1,2} = imgs{x, y, z+1};
               shifts{from,to} = c;
            end 
         end
      end
   end

   % calculate pairwise shifts

   for i = 1:np
      for j = 1:np
         if ~isempty(shifts{i,j})
            shifts{i,j} = fun(shifts{i,j}, param);
         end
      end
   end
   
   for i = 1:np
      for j = 1:np
         if ~isempty(shifts{i,j}) && isempty(shifts{j,i})
            shifts{j,i} = -shifts{i,j};
         end
      end
   end
   
   %var2char(shifts)

   % assemble shifts from primary directions

   start = getParameter(param, 'anchor', ceil(si/2));
   cstart = num2cell(start);

   sh = cell(si);
   sh{cstart{:}} = zeros(1,dim);

   sh = primaryShifts(sh, shifts, start, dim, si);
   sh = zero2first(sh);

end


function sh = primaryShifts(sh, shifts, start, dim, si)

   % loop over nearest neighbours for shifts
   rec = {};
   for d = 1:dim
      if start(d) < si(d) 
         neigb = start;
         neigb(d) = start(d) + 1;
         neigb = num2cell(neigb);
         if isempty(sh{neigb{:}})
            sh = findPrimaryShift(sh, shifts, [neigb{:}], dim, si);
            rec = {rec{:}, [neigb{:}]}; %#ok<CCAT>
         end
      end
      
      if start(d) > 1 
         neigb = start;
         neigb(d) = start(d) - 1;
         neigb = num2cell(neigb);   
         if isempty(sh{neigb{:}})
            sh = findPrimaryShift(sh, shifts, [neigb{:}], dim, si);
            rec = {rec{:}, [neigb{:}]}; %#ok<CCAT>
         end
      end
   end
   
   % loop over nearest neighbours for recursion
   for n = 1:length(rec)
      sh = primaryShifts(sh, shifts, rec{n}, dim, si);
   end
   
end


function sh = findPrimaryShift(sh, shifts, pos, dim, si)

   % find neighbours in different directions
   to = imsub2ind(si, pos);
   fr = cell(1,dim);

   for d = 1:dim
      if pos(d) < si(d)
         neigb = pos;
         neigb(d) = neigb(d) + 1;
         neigb = num2cell(neigb);
         if ~isempty(sh{neigb{:}})
            fr{d} = neigb;
         end
      end

      if isempty(fr{d}) && pos(d) > 1
         neigb = pos;
         neigb(d) = neigb(d) - 1;
         neigb = num2cell(neigb);
         if ~isempty(sh{neigb{:}})
            fr{d} = neigb;
         end
      end
   end
   
   % firs non-zero shift
   sfi = find(~cellfun(@isempty, fr), 1);
   if isempty(sfi)
      error('alignImagesOnGrid; inconsistent shifts!');
   end
   
   fri = fr{sfi};
   sfi = shifts{sub2ind(si, fri{:}), to} + sh{fri{:}};
   
   % calculate primary shift
   shi = zeros(1,dim);
   for d = 1:dim
      if ~isempty(fr{d})
         frd = fr{d};
         s = shifts{sub2ind(si, frd{:}), to} + sh{frd{:}};
         shi(d) = s(d);
      else
         shi(d) = sfi(d);
      end  
   end
   
   cpos = num2cell(pos);
   sh{cpos{:}} = shi;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% global shift optimization

function sh = alignByGlobal(imgs, fun, param)

   % determine pairs 

   si = size(imgs); dim = ndims(imgs{1});
   si = padright(si, 3, 1);

   np = si(1) * si(2) * (si(3)-1) + si(1) * (si(2)-1) * si(3) + (si(1)-1) * si(2) * si(3);
   shifts(np) = struct('from', -1, 'to', -1, 'shift', zeros(1,dim));

   np = 0;
   for x = 1:si(1)
      for y = 1:si(2)
         for z = 1:si(3)
            if x < si(1)
               np = np + 1;
               shifts(np).shift = {imgs{x,y,z}; imgs{x+1, y, z}}; 
               shifts(np).from = imsub2ind(si, [x,y,z]);
               shifts(np).to = imsub2ind(si, [x+1,y,z]);
            end
            if y < si(2)
               np = np + 1;
               shifts(np).shift = {imgs{x,y,z}, imgs{x, y+1, z}};
               shifts(np).from = imsub2ind(si, [x,y,z]);
               shifts(np).to = imsub2ind(si, [x,y+1,z]);
            end
            if z < si(3)
               np = np + 1;
               c{1,1,1} = imgs{x,y,z}; c{1,1,2} = imgs{x, y, z+1};
               shifts(np).shift = c;
               shifts(np).from = imsub2ind(si, [x,y,z]);
               shifts(np).to = imsub2ind(si, [x,y,z+1]);
            end 
         end
      end
   end
   

   % calculate pairwise shifts

   for i = 1:np
      shifts(i).shift = fun(shifts(i).shift, param);
   end
 
   % globally optimize pairwise shifts

   [~, ic] = optimizePairwiseShifts(shifts);
   ic = round(ic);
   
   % transform consistent shifts to shifts 
   
   sh = num2cell(ic',2);
   sh = zero2first(sh);
   sh = reshape(sh, si);
   
end



