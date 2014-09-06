function c = padcell(c, n, vals)
%
% c = padcell(c, n, vals)
%

if ~iscell(c)
   error('padcell: expects cell as first input!');
end

c = c(:)';
l = length(c);

if n > l
   vals = vals(:)';
   if ~iscell(vals)
      vals = num2cell(vals);
   end

   c = [c, repmat(vals, 1, ceil((n-l)/length(vals)))];
end

c = c(1:n);

end