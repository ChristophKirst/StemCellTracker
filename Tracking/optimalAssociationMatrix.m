function result = optimalAssociationMatrix(cost, A)
%
% result = optimalAssociationMatrix(A, cost)
%
% description:
%    Optimizes the association matrix A given the cost matrix cost.
%    It is assumed that creation/deletion of objects are the last
%    row and column in cost.
%
% input:
%    cost      cost matrix
%    A         optional initial association matrix
%
% output:
%    result    optimal association matrix given the cost matrix
%
% reference: 
%    A shortest augmenting path algorithm for dense and sparse linear assignment problems
%    Jonker, R, Volgenant, A, Computing 1987 
%
% todo: 
%    use c version of the above reference
%
% See also: matchObjects 


if nargin == 1 || isempty(A)
   A = initAssociationMatrix(cost);
end

%% check for consistency
if ~checkAssociationMatrix(A) 
    disp('optimalAssociationMatrix: association matrix failed consistency check');
    return
end

%% find matching
finished = 0;
while ~finished
    [A, finished] = doOneMove(A,cost);
end

%% check for consistency
if ~checkAssociationMatrix(A) % check A for consistency
    disp('optimalAssociationMatrix: association matrix failed consistency check');
    return
end

result = A;

end


%% Association Matrix Helpers

function ch = checkAssociationMatrix(A)
   %
   % check association matrix A for consistency
   %

   osize = size(A,1) - 1;
   nsize = size(A,2) - 1;

   ch = 1;
   s = sum(A(:,1:nsize),1);
   if find(s(1:nsize)~=1),
      disp('Inconsistent association matrix: Columns: ');
      find(s(1:nsize)~=1)
      ch = 0;
   end;
   s = sum(A(1:osize,:),2);
   if find(s(1:osize)~=1),
      disp('Inconsistent association matrix: Rows: ');
      find(s(1:osize)~=1)
      ch = 0;
   end
   
end


function A = initAssociationMatrix(cost)
%
% assume last row/col of C corresponds to dummy particle, others are real
%

   osize = size(cost,1) - 1;
   nsize = size(cost,2) - 1;

   A = zeros(size(cost));
   for i=1:osize,
       % sort costs of real particles
       [~,srtidx] = sort(cost(i,:));
       % append index of dummy particle
       iidx = 1;
       dumidx = find(srtidx==(nsize+1));
       % search for available particle of smallest cost or dummy
       % particle must not be taken and cost must be less than dummy
       while and(sum(A(:,srtidx(iidx)))~=0, iidx<dumidx), 
           iidx = iidx + 1;                               
       end
       A(i,srtidx(iidx)) = 1;
   end

   % set dummy particle for columns with no entry
   s = sum(A,1);
   A(osize+1,s < 1) = 1;
   % dummy always corresponds to dummy
   A(osize+1,nsize+1) = 1;

end



%% Optimization


function [A, finished]=doOneMove(A,C)
%
% replaces single link that maximizes the reduction in total cost
%

   osize = size(A,1) - 1;
   nsize = size(A,2) - 1;

   % find unmade links with finite cost
   todo = intersect(find(A(1:osize,1:nsize)==0),find(C(1:osize,1:nsize)<Inf));


   % determine induced changes and reduced cost cRed for each
   % candidate link insertion

   [iCand,jCand] = ind2sub([osize nsize],todo);
   k = size(iCand);

   cRed = zeros(k);
   xCand = zeros(k);
   yCand = zeros(k);

   for ic=1:k,

       xCand(ic) = find(A(iCand(ic),:)==1);
       yCand(ic) = find(A(:,jCand(ic))==1);

       cRed(ic) = C(iCand(ic),jCand(ic)) + C(yCand(ic),xCand(ic));
       cRed(ic) = cRed(ic)-C(iCand(ic),xCand(ic))-C(yCand(ic),jCand(ic));
   end
   
   % find minimum cost and corresponding action
   [minc,mini] = min(cRed);

   % if minimum is < 0, link addition is favorable
   if minc < -1e-10,
       A(iCand(mini),jCand(mini)) = 1;
       A(yCand(mini),jCand(mini)) = 0;
       A(iCand(mini),xCand(mini)) = 0;
       A(yCand(mini),xCand(mini)) = 1;
       finished = 0;
   else
       finished = 1;
   end

end


