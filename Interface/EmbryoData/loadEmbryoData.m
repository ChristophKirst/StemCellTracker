function frames = loadEmbryoData(dirname, param)
%
% frames = loadEmbryoData(dirname, param)
%
% description: 
%    loads the segmentated embryo data time frames in a directory
%    returns array of objects
%    with fields (r, time, size, intensity, id) for each time frame
%
% input:
%    dirname   directory containing the inout files
%    param     .load.min     - skip a file if number of objects is <= than this (2)
%              .load.change  - skip a file if abs(#objs prev - # objs current)/#obj prev 
%                              is larger than this ([] ignores this test)
%              .filter       - filter applied to each frame ([] for none)
%              .print.load   - print messages during loading files (0)
%
% output:
%    frames    array of frames
%
% See also: loadEmbryoDataFile

% initialize

if nargin < 2
   param = [];
end

load_min    = getParameter(param, {'load', 'min'}, 2);
load_change = getParameter(param, {'load', 'change'}, []);

filter = getParameter(param, {'filter'}, []);

print_load  = getParameter(param, {'print', 'load'}, 0);


% find data file names
fn = dir(fullfile(dirname, '*.csv'));
if ~isempty(fn)
    frames(1) = Frame();
else
   error('loadEmbryoData: no .csv files to open in %s\n', dir);
end


% extract the objects in each file
j = 1;
length_prev = Inf;

for i=1:length(fn)
    
   dat = loadEmbryoDataFile(fullfile(dirname,fn(i).name), param);
      
   if ~isempty(filter)
      dat = filter(dat);
   end 

    % ignore data if segmentation failed
    n = length(dat.objects);
    if n > load_min
      if ~isempty(load_change)
         if length_prev == Inf
            length_prev = n;
         end
         change = abs(n - length_prev)/n;
         if change < load_change
            frames(j) = dat; %#ok<AGROW>
            j = j+1;
            length_prev = n;
         else
            if print_load
               fprintf('loadEmbryoData: did not load data in %s\nchange in objects is to large %g > %g\n', fn(i).name, n, change, load_change);
            end
         end
      else
         frames(j) = dat; %#ok<AGROW>
         j = j+1;
      end
    else
      if print_load
         fprintf('loadEmbryoData: did not load data in %s\nnumber of objects to small %d <= %d\n', fn(i).name, n, load_min);
      end   
    end
        
end

if print_load
   fprintf('loadEmbryoData: loaded %d frames from %d files in %s\n\n', length(frames), length(fn), dirname)
end


end

