function [sub, isize] = imfrmtReformatSubIndex(sub, isize, infrmt, outfrmt)
%
% sub = imfrmtReformatSubIndex(sub, isize, infrmt, outfrmt)
%
% description: 
%     takes sub indices of image of size isize and format infrmt and transforms
%     to corresponding sub indices in outfrmt
%
% input:
%     sub        sub indices as e.g. returned by imind2sub
%     isize      image size
%     infrmt     (optional) format of input image (imfrmtFormat(img))
%     outfrmt    (optional) format of output image 
%
% output:
%     sub        reformatted indices
%
% note:
%     upper case letters are mathematically in postive axis, lower case mathematically inverted axis


if nargin < 3 || isempty(infrmt)
   infrmt = imfrmtFormatSize(isize);
end
infrmt = imfrmtFormat(infrmt);

if nargin < 4 || isempty(outfrmt)
   outfrmt = 'XYZCT';
end
outfrmt = imfrmtFormat(outfrmt);


if length(unique(lower(infrmt))) ~= length(infrmt)
   error('imfrmtReformatSubIndex: labels appear more than once in input format: %s', infrmt)
end
if length(unique(lower(outfrmt))) ~= length(outfrmt)
   error('imfrmtReformatSubIndex: labels appear more than once in output format: %s', outfrmt)
end
if length(infrmt) ~= size(sub,2)
   error('imfrmtReformatSubIndex: sub indices dont match informat')
end


outfrmtl = lower(outfrmt);
infrmtl = lower(infrmt);
nout = length(outfrmt);

%reverse coordinates
[id, pos] = ismember(infrmtl, outfrmtl);
pos = pos(id);
id = find(id);

for i= 1:length(id)
   if infrmt(id(i)) ~= outfrmt(pos(i))
      sub(:,id(i)) = isize(id(i)) - sub(:,id(i)) + 1;
   end
end

% add extra dimensions
for i = 1:length(outfrmtl)
   sf = strfind(infrmtl, outfrmtl(i));
   if isempty(sf)
      infrmtl(end+1) = outfrmtl(i);  %#ok<AGROW> % add to the end as this would introduce a new dimensions of size 1 in the original array
   end
end


% find dimensions to squeeze 
for i = 1:length(infrmtl)
   sf = strfind(outfrmtl, infrmtl(i));
   if isempty(sf) % remove dimensions
      outfrmtl = [outfrmtl, infrmtl(i)]; %#ok<AGROW> %% appending it effectively removes dimensions of size 1
   end
end
% now out_format and in_format have same size
%infrmt
%outfrmt

%find permutation
per = ones(1, length(infrmtl));
for i = 1:length(infrmtl)
   per(i) = strfind(infrmtl, outfrmtl(i));
end
%per

% permute indices
sub = sub(:, per(1:nout));

if nargout > 1
   isize = isize(per(1:nout));
end

end