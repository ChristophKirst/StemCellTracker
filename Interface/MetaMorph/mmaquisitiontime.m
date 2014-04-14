function t = mmaquisitiontime(filename, format)
%
% t = mmaquisitiontime(filename, format)
%
% description: 
%    returns acquisition-time-local property of image generated in metamorph software
%
% input:
%    filename   name of image file or folder
%    format     (optional) 'str'='string' or 'num', 'vec', 'seconds' ('string')
%
% output:
%    t          time as string or cell of strings
%
% See also: imfinfo

if nargin < 2
   format = 'string';
end

if isdir(filename)
   fns = dir(fullfile(filename, '*.TIF'));
   fns = cellfun(@(x) fullfile(filename, x), {fns.name},'UniformOutput', false);
elseif exist(filename, 'file') == 2
   fns = {filename};
else
   fns = dir(filename);
   fns = cellfun(@(x) fullfile(filename, x), {fns.name},'UniformOutput', false);
end

t = cell(1, length(fns));
for f = 1: length(fns)
   info = imfinfo(fns{f});
   
   if ~isfield(info, 'ImageDescription')
      error('bric_getAquisitionTime: image does not have ImageDescription meta data!');
   end

   for s = 1:length(info) % its possible to have a tif stack, then info is a struct array
      tinfo = regexp(info(s).ImageDescription, '.*<prop id="acquisition-time-local" type="time" value="(?<time>\d+ \d+:\d+:\d+.\d+)"/>.*', 'names');
      t{f}{s} = [tinfo.time(1:4) '-' tinfo.time(5:6) '-' tinfo.time(7:end)];
   end
   
   clear info;
end

fun = [];
switch format
   case {'num', 'days', 'datenum'}
      fun = @(x) datenum(x, 'yyyy-mm-dd HH:MM:SS.FFF');
   case {'vec', 'vector', 'datevec'}
      fun = @(x) datevec(x, 'yyyy-mm-dd HH:MM:SS.FFF');
   case 'seconds'
      fun = @(x) 24 * 3600 * datenum(x, 'yyyy-mm-dd HH:MM:SS.FFF');
   case {'milliseconds', 'ms'}
      fun = @(x) 1000 * 24 * 3600 * datenum(x, 'yyyy-mm-dd HH:MM:SS.FFF');  
end

if ~isempty(fun)
   for f = 1:length(t)
      t{f} = cellfun(fun, t{f}); %, 'UniformOutput', false);
   end
end

if ~isdir(filename)
   t = t{1};
end

