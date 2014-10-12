function isize = imfrmtPermuteSize(isize, infrmt, outfrmt)
%
% img = imfrmtPermuteSize(img, in_format, out_format)
%
% description: 
%     takes size isize of an image with format infrmt and permutes it to the outfrmt
%     procudes error if non-songelton dimensions are lost 
%
% input:
%     isize      input image size
%     infrmt     (optional) format of input image (imfrmtFormat(img))
%     outfrmt    (optional) format of output image 
% output:
%     isize      reformatted image size
%

if nargin < 2 || isempty(infrmt)
   infrmt = imfrmtFormatSize(isize);
end
infrmt = imfrmtFormat(infrmt);

if nargin < 3 || isempty(outfrmt)
   outfrmt = 'XYZCT';
end
outfrmt = imfrmtFormat(outfrmt);


if length(unique(lower(infrmt))) ~= length(infrmt)
   error('imfrmtPermuteSize: labels appear more than once in input format: %s', infrmt)
end
if length(unique(lower(outfrmt))) ~= length(outfrmt)
   error('imfrmtPermuteSize: labels appear more than once in output format: %s', outfrmt)
end

% size is independent of orientation -> lower case
infrmt = lower(infrmt);
outfrmt = lower(outfrmt);

% add extra dimensions
for i = 1:length(outfrmt)
   sf = strfind(infrmt, outfrmt(i));
   if isempty(sf)
      infrmt(end+1) = outfrmt(i);  %#ok<AGROW> % add to the end as this would introduce a new dimensions of size 1 in the original array
   end
end


% find dimensions to squeeze 
for i = 1:length(infrmt)
   sf = strfind(outfrmt, infrmt(i));
   if isempty(sf) % remove dimensions only if it is of size 1 (squeeze)
      if isize(i) == 1
         outfrmt = [outfrmt, infrmt(i)]; %#ok<AGROW> %% appending it effectively removes dimensions of size 1
      else
         error('imfrmtReformatSize: non-trivial dimension %s not specified in output format %s!', infrmt(i), outfrmt);
      end
   end
end
% now out_format and in_format have same size

%find permutation
per = ones(1, length(infrmt));
for i = 1:length(infrmt)
   per(i) = strfind(infrmt, outfrmt(i));
end

%append ones to size
isize = padright(isize, length(outfrmt), 1);

% premute
isize = isize(per);

end
   