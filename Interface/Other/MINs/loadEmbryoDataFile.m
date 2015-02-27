function frame = loadEmbryoDataFile(fn, param)
%
% frame = loadEmbryoDataFile(fn, param)
%
% description: 
%    loads the segmentated embryo data from a single time frame
%    and stores the data in an array of Object classes
%
% input:
%    fn     filename of file contaning a comma separated table with fields
%           'X', 'Y', 'Z', 'Size', 'CH1-Avg', 'Cell ID'
%    param  .print.load   - print messages while loading the file (0)        
%
% output:
%    frame  Frame classes representing the data in fn  
%
% See also: loadEmbryoData

% send messages / info to user on first read 
persistent size_flag time_flag;


if nargin < 2
   param = [];
end

print_load  = getParameter(param, {'print', 'load'}, 0);



% labels for the relevant data 
% order should be: {CellID, X, Y, Z, (Size), (Intensity), ...} 
% size and intensity are optional, further entries can be used for
% typing of cells (i.e. Embryo Id)

data_label = {'Cell ID', 'X', 'Y', 'Z', 'Size', 'CH1-Avg'};


% time
[idx1, idx2] = regexpi(fn, 'frame=[0-9]+' ); 
time = str2double( fn(6+idx1:idx2) );

if print_load && isempty(time_flag)
    fprintf('loadEmbryoDataFile: inferred time or frame number = %d\nfrom file= %s\n', time, fn);
    time_flag = time;
end



% import data, check and find labels
data = importdata(fn,',', 1);

s = size(data.data,2);
if s < 8
    error('loadEmbryoDataFile: number of data columns is too small (<8): %d\nin file: %s\n', s, fn);  
elseif print_load && isempty(size_flag) && s < 12
   fpirntf('loadEmbryoDataFile: warning: number of data columns is not 12 but %d\nin file: %s\n', s, fn);
   size_flag = s;
end

labels = strsplit(data.textdata{1}, ',');
labels = strtrim(labels);
%if labels{end} == ''
%   labels = labels{1:end-1};
%end

[~, data_inds] = ismember(data_label,labels);

if ismember(0, data_inds)
   id = find(data_inds==0, 1);
   error('loadEmbryoDataFile: cannot find all labels %s\nin file: %s\n', char(data_label(id)), fn);
end

if length(data_inds) < 4
   error('loadEmbryoDataFile: number of data labels to small\nneed at least {id, x, y, z}\n in file: %s\n', fn);
end
  

% extract relevant data: (cellid, x, y, z, size, intensity)

data = data.data;
data = data(:,data_inds);

% convert to Object class

for n = size(data,1):-1:1 
   dd = data(n,:);
   objs(n) = Object('id', dd(1), 'r', dd(2:4)', 'volume', dd(5), 'intensity', dd(6), 'time', time);
end

frame = Frame('t', time, 'objects', objs);

end
