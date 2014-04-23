function tags =findtags(tfrmt)

fns = tagformat2files(tfrmt);
tags = name2tags(tfrmt, fns);

% tr = cell(size(vals,1),1);
% for i = 1:size(vals,1);
%    switch types{i}
%       case 's'
%          tr{i} = unique(vals(i,:));
%       otherwise
%          tr{i} = unique([vals{i,:}]);
%    end
% end

end