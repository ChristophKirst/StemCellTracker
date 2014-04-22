function [objects, varargout] = label2DataObjects(imglab, img, stats, statnames, param)
%
% objects = label2DataObjects(imglab, img, stats, param)
% [objects, stats] = label2DataObjects(...)
%
% description:
%    converts a labeld image and associated data to an array of DataObject classes
%    used for tracking and further analysis
%
% input:
%    imglab    labeled image (2D/3D)
%    img       original grayscale image to measure basic properties from (2D/3D)
%    stats     (optional) previously calcualted statistics
%    statnames (optional) names of measures to actually save in the DataObject ('existing' = all names in the stats struct, 'all' = calculate all possible statistics)
%    param     (optional) parameter struct with entries:
%              .time      time for objects (0)
%              .rescale   rescale coordinates r by this factor ([1, 1(, 1)])
%              .method    how to calcualte the intensity in Object tracking, a string of any function, 'none' = dont calcualte ('median')
%  
% output:
%    objects array of DataObjects each representing one of the labels
%    stats   (optional) updated statistics
%
% See also: Object, DataObject, imstatistics

if nargin < 3
   stats = struc();
end

if nargin < 4
   statnames = 'existing';
   param = [];
else
   if isstruct(statnames)
      param = statnames;
      statnames = 'existing';
   else
      if nargin < 5
         param = [];
      end
   end
end

if isempty(statnames)
   statnames = 'existing';
end

if ~iscell(statnames)
   statnames = {statnames};
end

if any(ismember(statnames, 'existing'))
   statnames = fieldnames(stats);
end

%fprintf('label2DataObjects: calculating statistiics...\n');
statnames = {statnames{:}}; %#ok<CCAT1>
statnames = [statnames, {'Volume', 'Centroid', 'PixelIdxList'}]; % we definately need these for the tracking data!
stats = imstatistics(imglab, stats, statnames, img);


time = getParameter(param, 'time', 0);
if isempty(time)
   time = 0;
end

dim = ndims(imglab);

rescale = getParameter(param, 'rescale', 1);
if isempty(rescale)
   rescale = 1;
end
rescale = rescale(:)';
rescale = padright(rescale, dim, 'circular');

meth = getParameter(param, 'method', 'median');

nobj = length(stats);


% Start creating DataObjects
objects(nobj) = DataObject();

%id
ids = {stats.PixelIdxList};
ids = cellfun(@first, ids);
ids = imglab(ids);
ids = num2cell(ids);

% length(ids)
% length(objects)
% 
% size(ids)
% size(objects)

[objects.id] = ids{:};

%time
[objects.time] = deal(time);

%r
rr = repmat(rescale', 1, nobj) .* [stats.Centroid];
rr = num2cell(rr,1);
[objects.r] = rr{:};

%volume
[objects.volume] = stats.Volume;

if ~strcmp(meth, 'none')
   fun = str2func(meth);
 
   val = cellfun(@(x) fun(img(x)), {stats.PixelIdxList}, 'UniformOutput', false);
   [objects.intensity] = val{:};
   [stats.([upper(meth(1)), meth(2:end), 'Intensity'])] = val{:};
end

%data entries
statsnew = rmfield(stats,  {'Volume', 'Centroid', 'PixelIdxList'});
statsnew = num2cell(statsnew);
[objects.data] = statsnew{:};


%compressed segmentation
segs(nobj) = Segment();
isize = size(imglab);

for i = 1:nobj
   %fprintf('label2DataObjects: compressing label: %g / %g\n', i, nobj);
   segs(i).store(stats(i).PixelIdxList, isize);
end

segs = num2cell(segs);
[objects.segment] = segs{:};


% return updated statistics
if nargout > 1
   varargout{1} = stats;
end

end


