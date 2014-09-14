function [regs, ids] = stitchImagesOverlapRegions(ashifts, isizes)
%
% [regs, ids] = stitchImagesOverlapRegions(ashifts, isizes)
%
% description:
%    splits a set of co-axial images aligned by shifts ashifts and sizes isizes into a minimal set of
%    non-overlaping regions used by stiching routines
%
% input:
%    ashifts   cell array of shifts as row vectors
%    isizes    cell array of image sizes as row vectors
%
% output:
%    regs      cell array of non-overlapping regoins as [pos1; pos2] where pos1, pos2 are row vectors of corners with pos1(i) < pos2(i)
%    ids       cell array of ids of images in each regoin
%
% See also: stitchImages

n = numel(ashifts);

regs = {};
ids  = {};
for i = 1:n
   [regs, ids] = addOverlapRegion(regs, ids, ashifts{i}, isizes{i}, i); 
   
%    disp(i)
%    var2char({'regs', regs})
%    var2char({'ids',  ids })
%    disp ====================================
   
   
end

end
   

function [regsnew, idsnew] = addOverlapRegion(regs, ids, shift, isize, id)
   
   dim = length(isize);

   % try too add the full new region first
   regsadd = {[shift + 1; shift + isize]};

   % these are non-overlapping rectangles to bechecked for overlap with regsadd
   regscheck = regs;
   idscheck  = ids;

   % rectables that will now have any further overlap with regsadd
   regsnew   = {};
   idsnew    = {};

   i = 1;
   while ~isempty(regscheck) && ~isempty(regsadd)
 
      % chek the next one in list
      rc = regscheck{1};
      
      %check if overlap with any of the regoins to be added
      found = false;
      for a = 1:length(regsadd)
         
         ra = regsadd{a};
         ov = findOverlap(rc, ra);
         
         if ~isempty(ov) % overlapping region is save to add to new list and remove form check and add list after splitting
            
            
            % split region to check
            split = splitRegion(rc, ov, dim);
            
%             var2char({'regsadd', regsadd})
%             var2char({'regscheck', regscheck})
%             var2char({'ra', ra, 'rc', rc})
%             var2char({'split',   split})

            % add id of immage to add to overlapping region -> cannot overlap with other regoin -> safe to put into new list
            regsnew{i} = split{1}; %#ok<AGROW>
            idsnew{i}  = [idscheck{1}, id]; %#ok<AGROW>
            i = i + 1;
            
            % other non-overlapping regoins from split need to beckeched anew
            regscheck = [split(2:end), regscheck(2:end)];
            idscheck  = [repmat(idscheck(1), 1, length(split)-1), idscheck(2:end)];
            
            
            % split added region
            split = splitRegion(ra, ov, dim);
          
            % add all non-verlapping regions to the tob be added list
%             var2char({'regsadd before: ' , regsadd})
%             var2char({'split ra: ' , split})
             regsadd = [regsadd(1:a-1), split(2:end), regsadd(a+1:end)];
%             var2char({'regsadd after: ' , regsadd})
            
            % start cehcking anew
            found = true;
            break;
            
         end
      end
      
      if ~found   % -> check region not overlapping -> add to new list and remove from checklist
         regsnew{i} = rc; %#ok<AGROW>
         idsnew{i}  = idscheck{1}; %#ok<AGROW>
         i = i + 1;
         
         regscheck = regscheck(2:end);
         idscheck  = idscheck(2:end);
      end
      

%       var2char({'check', regscheck})
%       var2char({'add', regsadd})
%       var2char({'new', regsnew})
%       disp -----------------------------------------
    
   end

   regsnew = [regsnew, regscheck, regsadd];
   idsnew  = [idsnew,  idscheck,  num2cell(id * ones(1,length(regsadd)))];


end



%splits a rectangle r given the overlap region o \subset r
function split = splitRegion(r, o, dim)
   if dim ==2
      split = splitRegoin2D(r,o);
   else
      split = splitRegoin3D(r,o);
   end
end

function split = splitRegoin2D(r, o)
   r1 = r(1,:); r2 = r(2,:);
   o1 = o(1,:); o2 = o(2,:);
         
   split = {o};
   
   % split x-axis
   if r1(1) < o1(1)
      split = [split, {[r1; o1(1)-1, r2(2)]}];
   end
   
   if o2(1) < r2(1)
      split = [split, {[o2(1)+1, r1(2); r2]}];
   end
   
   % split y-axis
   if r1(2) < o1(2)
      split = [split, {[o1(1), r1(2); o2(1), o1(2)-1]}];
   end
   
   if o2(2) < r2(2)
      split = [split, {[o1(1), o2(2)+1; o2(1), r2(2)]}];
   end
end

function split = splitRegoin3D(r, o)
   r1 = r(1,:); r2 = r(2,:);
   o1 = o(1,:); o2 = o(2,:);
   
   split = {o};
   
   % split x-axis
   if r1(1) < o1(1)
      split = [split, {[r1; o1(1)-1, r2(2), r2(3)]}];
   end
   
   if o2(1) < r2(1)
      split = [split, {[o2(1)+1, r1(2), r1(3); r2]}];
   end
   
   % split y-axis
   if r1(2) < o1(2)
      split = [split, {[o1(1), r1(2), r1(3); o2(1), o1(2)-1, r2(3)]}];
   end
   
   if o2(2) < r2(2)
      split = [split, {[o1(1), o2(2)+1, r1(3); o2(1), r2(2), r2(3)]}];
   end
    
   % split z-axis
   if r1(3) < o1(3)
      split = [split, {[o1(1), o1(2), r1(3); o2(1), o2(2), o1(3)-1]}];
   end
   
   if o2(3) < r2(3)
      split = [split, {[o1(1), o2(2), o2(3)+1; o2(1), o2(2), r2(3)]}];
   end
   
end


function ov = findOverlap(a,b)
   a1 = a(1,:); a2 = a(2,:);
   b1 = b(1,:); b2 = b(2,:);
   
   ov = zeros(2,length(a1));
   ov(1,:) = max(a1,b1);
   ov(2,:) = min(a2,b2);

   if any(ov(2,:)-ov(1,:) < 0)
      ov = [];
   end
   
   
  % {var2char(a), var2char(b), var2char(ov)}

end


