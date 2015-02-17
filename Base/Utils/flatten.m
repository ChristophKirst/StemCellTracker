function C = flatten(A)
% 
% C = flatten(A)
%
% description: 
%     flattens cell array to top level

C = {};
for i=1:numel(A)  
    if(~iscell(A{i}))
       C = [C,A{i}]; %#ok<AGROW>
    else
       Ctemp = flatten(A{i});
       C = [C,Ctemp{:}]; %#ok<AGROW>
    end
end