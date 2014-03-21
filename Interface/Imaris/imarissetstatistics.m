function istat = imarissetstatistics(varargin)
%
% istat = imarissetstatistics(name, values)
% istat = imarissetstatistics(objectname, ...)
% istat = imarissetstatistics(object, ...)
% istat = imarissetstatistics(imaris, ...)
%
% description:
%    add statistic values for a surfaces object in Imaris, if struct is
%    passed as values all scalar entries are added
%
% input:
%    name           name of the statistics, or cell of names
%    values         array or cell of values for each surface, or struct with statistics
%
%    objectname     name of Imaris surface
%    object         Imarise ISurface surface
%    imaris         Imaris application instance
%
% output:
%    istat   statistics
%
% note: 
%    standard Factors used by Imaris and here:
%        'Category'   = 'Surface'
%        'Channel'    = ''
%        'Collection' = ''
%        'Time'       = '1'
%
% See also: imarisset

[imaris, varargin, nargin] = imarisvarargin(varargin);

if nargin < 1 || nargin > 3
   error('imarissetstatistics: expect 1-3 input arguments');
end

%add = 0;
if ischar(varargin{1}) && nargin > 1 && ischar(varargin{2})
   surfacename = varargin{1};
   isurface = imarisgetobject(imaris, surfacename, 'Surfaces');
   if isempty(isurface)
      isurface = imariscreatesurface(imaris, surfacename);
      %add = 1;
   end
   varargin = varargin(2:end);
   %nargin = length(varargin);
elseif isimaristype(imaris, varargin{1}, 'Surfaces')
   isurface = varargin{1};
   varargin = varargin(2:end);
   %nargin = length(varargin);
else % try to get selected surface
   isurface = imarisgetcurrentobject('Surfaces');
   if isempty(surface)
      error('imarissetstatistics: select a valid surfaces object!')
   end
end

if isempty(surface)
   error('imarissetstatistics: cannot determine surfaces object!')
end

try
   nsurfaces = isurface.GetNumberOfSurfaces();
   istat = isurface.GetStatistics();
catch %#ok<CTCH>
   error('imarissetstatistics: error in retrieving current statistics!')
end

if nsurfaces == 0
   error('imarissetstatistics: error no surfaces to set statistics!')
end


if nargin == 1 && isstruct(varargin{1}) % statistics struct
   stats = varargin{1};
   if length(stats) ~= nsurfaces
      error('imarissetstatistics: inconsistent size of statisticas and number of surfaces!')
   end
   
   allnames = fieldnames(stats);
   names = {};
   values = {};
   for n = 1:length(allnames)
      if isscalar(stats(1).(allnames{n}))
         names = [names, allnames(n)]; %#ok<AGROW>
         values = [values, {stats.(allnames{n})}]; %#ok<AGROW>
      end
   end
   
elseif nargin == 2 % names value paris
   names = varargin{1};
   if ischar(names)
      names = {names};
   end
   if ~iscellstr(names)
      error('imarissetstatistics: expects name or cell array of names!')
   end
   
   values = varargin{2};
   if ~iscell(values)
      values = {values(:)};
   end
end


nstats = length(names);
if nstats ~= length(values)
   error('imarissetstatistics: nunmber of names does not match number of statitics!')
end

for i = 1:nstats
   if ~isvector(values{i}) || length(values{i}) ~= nsurfaces
      error('imarissetstatistics: dimension mistmatch!')
   end
end
      

%  factors
vFactorNames = {'Category', 'Channel', 'Collection', 'Time'}';
vFactors = {'Surface', '', '', '1'}';

 for s = 1:nstats
    statNames   = cell(nsurfaces, 1);
    statUnits   = cell(nsurfaces, 1);
    statFactors = cell(size(vFactors, 1), nsurfaces);
    statIds     = 0:(nsurfaces-1);
    %statIds     = 1:nsurfaces;
    for j = 1:nsurfaces
       statNames{j}      = names{s};
       statUnits{j}      = 'um';
       statFactors(:, j) = vFactors;
       %statIds(j)        = j;
    end
    statValues  = values{s};
    
    try
       isurface.AddStatistics(statNames, statValues, statUnits, statFactors, vFactorNames, statIds');
    catch %#ok<CTCH>
       error('imsetstatistics: error while setting the statistics!')
    end
 end














