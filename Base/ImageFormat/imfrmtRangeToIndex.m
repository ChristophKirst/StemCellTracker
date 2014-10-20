function idx = imfrmtRangeToIndex(isize, ifrmt, varargin)
%
% trg = imfrmtRangeToIndex(isize, ifrmt, rgs)
%
% description: 
%     converts range specs to array of indices given the size and format
%
% input:
%     isize      the size of the data array
%     ifrmt      the reference format
%     rgs        coordinate ranges
%
% output:
%     id         array of indices as array of indices corresponding to the ranges
%
% note: array is of shape corresponding to the format
% 
% See also: imfrmtIndexToCoordinate


% generate indices
rgs = parseParameter(varargin);

tnames = fieldnames(rgs);

nfrmt = length(ifrmt);

args = cell(1,nfrmt);
for i = 1:nfrmt
   id = find(strcmpi(tnames, ifrmt(i)), 1);
   if ~isempty(id) 
      % get restricted range
      v = rgs.(tnames{id});
      if iscell(v)
         v = cell2mat(v);
      end
   else  % get full range
      v = 1:isize(i);
   end  
   
   % flip if required
   if ~isempty(id) && ifrmt(i) ~= tnames{id}
      v = isize(i) - v +1;
   end
   
   args{i} = v;
end

gr = ndgridc(args{:});
si = size(gr{1});

gr = cellfunc(@(x) x(:), gr);
gr = cell2mat(gr);
% 
% isize
% ifrmt
% gr

idx = imsub2ind(isize, gr);
idx = reshape(idx, si);

end








   