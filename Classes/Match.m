classdef Match < handle
%
% Match is class storing data from matching two frames
%
   properties      
      n0 = 0;     % number of pre objects 
      n1 = 0;     % number of post objects
      match = []; % links of the form i -> match(i) =-1 if object is deleted
      
      objects0 = []; % pre objects 
      objects1 = []; % post objects
      
      cost = [];  % optional cost matrix used for the linking
   end

   methods
      function obj = Match(match, objects0, objects1, cost)
      %
      % obj = Match(match, n1, cost)
      % obj = Match(match, objects0, objects1, cost)
      %
      % input: 
      %   match     links of the form i -> match(i) =-1 if object is deleted
      %   n1        number of objects in frame 1 (= length(match) if not supplied)
      %   objects*  arrays of objects to which match refers to
      %   cost      optional cost matrix
      %
      
         if nargin > 0
            obj.n0 = length(match);
            obj.match = match;
      
            if isnumeric(objects0)

               obj.n1 = objects0;

               if nargin > 2
                  obj.cost = objects1;
               end

            else

               obj.objects0 = objects0;   
               if obj.n0 ~= length(objects0)
                  error('number of matching ids n0=%d and objects0=%d inconsistent', obj.n0. length(objects0))
               end       

               obj.objects1 = objects1;
               obj.n1 = length(objects1);

               if nargin > 3
                  obj.cost = cost;
               end

            end
         end
      end 
      
      function d = dim(obj)
      %
      % d = dim(obj)
      %
      % description:
      % returns spatial dimnsion of the objects
      %
         d = obj.objects0(1).dim;
         
      end
      
       
      function pairs = toArray(obj, creation_id, deletion_id)
      %
      % pairs = toArray(obj, creation_id, delete_id)
      %
      % description:
      %    returns array of paried indices corresponding to the matches
      %    [i match(i)] taking into account 
      %    creation of new objects as [creation_id j]
      %    deletion of objects as [i deletion_id]
      %
      % input:
      %    creation_id  number to indicate creation of objects ([] = n1 + 1)
      %    deletion_id  number to indicate deletion of objects ([] = n0 + 1)
      %
      % output:
      %    pairs    array of pairs [i match(i)]
      
         % initialize 
  
         if nargin < 2
            creation_id = -1;
         end
         if nargin < 3
            deletion_id = -1;
         end
    
         if isempty(creation_id)
            creation_id = obj.n1+1;
         end
         if isempty(deletion_id)
            deletion_id = obj.n0+1;
         end

         % calcualte pairs

         i = 1:obj.n0;
         j = obj.match'; 
         j(obj.match<0) = deletion_id;

         j = [j setdiff(1:obj.n1, j)];
         i = [i creation_id * ones(1,length(j)-length(i))];

         pairs = [i' j'];

      end
      

      function [X0, X1] = toCoordinates(obj)
      %
      % [X0, X1] = toCoordinates(obj)
      %
      % description:
      %    returns the coordinates of the objects of the matched pairs
      %
      % output:
      %    X0,X1    matched points as columns
      %
      
         if isempty(obj.objects0) || isempty(obj.objects1)
            error('cannot calculate positions without objects!')
         end
      
         % get pairs
         pre = find(obj.match > 0);
         post = obj.match(pre);
         

         X0 = [obj.objects0(pre).r];
         X1 = [obj.objects1(post).r];

      end    

      
      function c = matchingCost(obj, cost)
      %
      % c = matchingCost(obj, cost)
      %
      % description:
      %    compute cost for the match
      %
      % input:
      %    cost    cost matrix if not stored in Match object
      %

         cmat = obj.cost;
         if isempty(cmat) && nargin > 1
            cmat = cost;
         end
         if isempty(cmat)
            error('cannot calculate cost without cost matrix!');
         end

         c = 0;
         for i = 1:length(obj.match)
            j = obj.match(i);
            if j < 0
               j = obj.n0 + 1;
            end
            c = c + cost(i, j);
         end
         % add creation cost
         indx = setdiff(1:obj.n1, obj.match);
         c = c + sum(cost(nsize+1, indx));

      end
         
  
            
      function stats = statistics(obj, cost)
      %
      % stats = statistics(obj, data)
      %
      % description:
      %    returns some statistics of the trajectory
      %
      % input:
      %    cost      cost matrix
      %
      
         if nargin < 2
            stats = statisticsMatch(obj);
         else
            stats = statisticsMatch(obj, cost);
         end
 
      end
      

   end % methods
end % classdef
      
      
      
      
      