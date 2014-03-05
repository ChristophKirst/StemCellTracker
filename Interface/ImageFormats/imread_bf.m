function [data, metatdata] = imread_bf(name, param) 
%
% [data, metatdata] = imread_bf(name, param) 
%
% description:
%    import bio-format file by opening the loci importer gui
%
% input:
%    name   (optional) filename
%    param  (optional) struct with optional entries (if not present gui is opened, [] imports everything)
%           .series    ids of series to import (if array data is cell array) ([])
%           .time      ids of time frames to import as array [tid1, tid2, ...] ([] = all)
%           .channel   ids of channels to import as array [cid1, cid2, ...] ([] = all)
%           .z /.l     pixel range [min, max] in z / l direction ([])
%           .x /.h     pixel range [min, max] in x / h direction ([])
%           .y /.w     pixel range [min, max] in y / w direction ([])
%           .metadata  return also meta data if true
%           .gui       open gui to select series if true (full image series is improted)
%
% output:
%    data     image data
%    metadata (optional) image meta data 
%
% note: the loci_tools.jar must be present in the javaclasspath
%       try bfinitialize or ijintialize to set this path
%
% See also: bfinitialize, imread, ijimage2mat


% Prompt for a file if not input
if nargin == 0
  [file, path] = uigetfile(bffileextensions, 'Choose a file to open');
  name = [path file];
end
if ~exist(name, 'file')
   error('imread_bf: file does not exists: %s', name)
end

if nargin < 2
   param.gui = true;
end

gui = getParameter(param, {'gui'}, 0);

if gui
   data = importlocigui(name);
else
   [data, metatdata] = importlociseries(name, param);
end

end



% open file using java gui interface   
function data = importlocigui(name)

%get a loci importer
lociimp = loci.plugins.LociImporter();
imp = loci.plugins.in.Importer(lociimp);

options = imp.parseOptions(name);
if (lociimp.canceled) 
   return
end
   
process = loci.plugins.in.ImportProcess(options);
imp.showDialogs(process);

if (lociimp.canceled)
   return
end

displayHandler = loci.plugins.in.DisplayHandler(process);
displayHandler.displayOriginalMetadata();
displayHandler.displayOMEXML();

reader = loci.plugins.in.ImagePlusReader(process);
imps = imp.readPixels(reader, options, displayHandler);

%displayHandler.displayImages(imps);
%displayHandler.displayROIs(imps);

imp.finish(process);

% convert image to data
if length(imps) ==1 
   data  = ijimage2mat(imps(1));
else
   for i = length(imps):-1:1
      data{i} = ijimage2mat(imps(1));
   end
end

end







function [data, metadata] = importlociseries(name, param)

% initialize logging
loci.common.DebugTools.enableLogging('INFO');

% Get the channel filler
r = loci.formats.ChannelFiller();
r = loci.formats.ChannelSeparator(r);
%if stitchFiles
%    r = loci.formats.FileStitcher(r);
%end
r.setMetadataStore(loci.formats.MetadataTools.createOMEXMLMetadata());
r.setId(name);

numSeries = r.getSeriesCount();

s = getParameter(param, {'series'}, []);

if isempty(s)
   s = 1:numSeries;
end

if ~iscell(s) && numel(s) > 1
   s = mat2cell(s(:), ones(1, numel(s)));
end

if iscell(s)
   par = param;
   for i = length(s):-1:1
      par.series = i;
      [d, m] = importlociseries(name, par);
      data{i} = d; metadata{i} = m;
   end
   return
end

if s < 1 || s > numSeries
   error('imread_bf: series index out of range!');
end

%%% read a single series 

r.setSeries(s - 1);
numImages = r.getImageCount();

% dimensions
sizeZ = r.getSizeZ();
sizeC = r.getSizeC();
sizeT = r.getSizeT();
sizeX = r.getSizeX();
sizeY = r.getSizeY();

tids = getParameter(param, {'time'}, []);
if isempty(tids)
   tids = 1:sizeT;
else
   tids = tids(:)';
   if max(tids) > sizeT || min(tids) < 1
      error('imread_bf: time indices out of range !')
   end
   sizeT = length(tids);
end

cids = getParameter(param, {'channel'}, []);
if isempty(cids)
   cids = 1:sizeT;
else
   cids = cids(:)';
   if max(cids) > sizeC || min(cids) < 1
      error('imread_bf: channel indices out of range !')
   end
   sizeC = length(cids);
end

if isentry(param, 'z')
   zids = getParameter(param, {'z'}, []);
else
   zids = getParameter(param, {'l'}, []);
end
if isempty(zids)
   zids = 1:sizeZ;
else
   zids = zids(:)';
   if max(zids) > sizeT || min(zids) < 1
      error('imread_bf: z indices out of range !')
   end
   sizeZ = length(zids);
end

if isentry(param, 'x')
   xrange = getParameter(param, {'x'}, []);
else
   xrange = getParameter(param, {'h'}, []);
end
if isempty(xrange)
   x = 1; h = sizeX;
elseif length(xrange) == 1
   x = xrange(1); h = sizeX - xrange(1) + 1;
elseif length(xrange) ==2
   x = xrange(1); h = xrange(2) - xrange(1) + 1;
else
   error('imread_bf: x/h should be of the form min or [min, max]');
end
if x < 1 || x + h > sizeX + 1
   error('imread_bf: x/h [min, max] out of bounds');
else
   sizeX = h;
end
   

if isentry(param, 'y')
   yrange = getParameter(param, {'y'}, []);
else
   yrange = getParameter(param, {'w'}, []);
end
if isempty(yrange)
   y = 1; w = sizeX;
elseif length(yrange) == 1
   y = yrange(1); w = sizeX - yrange(1) + 1;
elseif length(xrange) ==2
   y = yrange(1); w = yrange(2) - yrange(1) + 1;
else
   error('imread_bf: x/h should be of the form min or [min, max]');
end
if y < 1 || y + w > sizeY + 1
   error('imread_bf: y/w [min, max] out of bounds');
else
   sizeY = w;
end



% Test plane size

planeSize = loci.formats.FormatTools.getPlaneSize(r, h, w);
%planeSize = loci.formats.FormatTools.getPlaneSize(r);
%end

if planeSize/(1024)^3 >= 2,
    error('imread_bf: image larger than 2GB! try opening i');
end

if 1 < s &&  s > numSeries
   error('imread_bf: series index out of bound. max series = %d', numSeries);
end


%pixelType = r.getPixelType();
%bpp = loci.formats.FormatTools.getBytesPerPixel(pixelType);
%bppMax = power(2, bpp * 8);
%imageList = cell(numImages, 2);
%colorMaps = cell(numImages);


data = zeros([sizeX, sizeY, sizeZ, sizeC, sizeT]);


fprintf('Reading series #%d', s);
for i = 1:numImages
   if mod(i, 72) == 1
      fprintf('\n    ');
   end
   fprintf('.');
   %arr = bfGetPlane(r, i, varargin{:});

   % retrieve color map data
%    if bpp == 1
%       colorMaps{s, i} = r.get8BitLookupTable()';
%    else
%       colorMaps{s, i} = r.get16BitLookupTable()';
%    end
%    
%    warning off %#ok<WNOFF>
%    if ~isempty(colorMaps{s, i})
%       newMap = single(colorMaps{s, i});
%       newMap(newMap < 0) = newMap(newMap < 0) + bppMax;
%       colorMaps{s, i} = newMap / (bppMax - 1);
%    end
%    warning on %#ok<WNON>
   
   
   % build an informative title for our figure
%    label = name;
%    if numSeries > 1
%       seriesName = char(r.getMetadataStore().getImageName(s - 1));
%       if ~isempty(seriesName)
%          label = [label, '; ', seriesName];
%       else
%          qs = int2str(s);
%          label = [label, '; series ', qs, '/', int2str(numSeries)];
%       end
%    end
%    if numImages > 1
%       qi = int2str(i);
%       label = [label, '; plane ', qi, '/', int2str(numImages)];
%       if r.isOrderCertain()
%          lz = 'Z';
%          lc = 'C';
%          lt = 'T';
%       else
%          lz = 'Z?';
%          lc = 'C?';
%          lt = 'T?';
%       end
%       zct = r.getZCTCoords(i - 1);
%       sizeZ = r.getSizeZ();
%       if sizeZ > 1
%          qz = int2str(zct(1) + 1);
%          label = [label, '; ', lz, '=', qz, '/', int2str(sizeZ)];
%       end
%       sizeC = r.getSizeC();
%       if sizeC > 1metadata = r.getSeriesMetadata();
%          qc = int2str(zct(2) + 1);
%          label = [label, '; ', lc, '=', qc, '/', int2str(sizeC)];
%       end
%       sizeT = r.getSizeT();
%       if sizeT > 1
%          qt = int2str(zct(3) + 1);
%          label = [label, '; ', lt, '=', qt, '/', int2str(sizeT)];
%       end
%    end
   
   % save image plane and label into the list
   %imageList{i, 1} = arr;
   %imageList{i, 2} = label;
      
   zct = r.getZCTCoords(i-1) + 1;
   z = find(zids == zct(1));
   c = find(cids == zct(2));
   t = find(tids == zct(3));
   
   if ~isempty(z) && ~isempty(c) && ~isempty(t)
      %fprintf('%g %g %g\n', z, c, t);
      data(:,:,z,c,t) = bfgetplane(r, i, x, y, h, w);
   end
   
end
fprintf('\n');


% save images and metadata into our master series list
%result{s, 1} = imageList;

% extract metadata table for this series
%result{s, 2} = r.getSeriesMetadata();
%result{s, 3} = colorMaps;
%result{s, 4} = r.getMetadataStore();

metadata = r.getSeriesMetadata();

r.close();

end


function I = bfgetplane(r, iPlane, x, y, h, w)
% Get pixel type
pixelType = r.getPixelType();
bpp = loci.formats.FormatTools.getBytesPerPixel(pixelType);
fp = loci.formats.FormatTools.isFloatingPoint(pixelType);
sgn = loci.formats.FormatTools.isSigned(pixelType);
little = r.isLittleEndian();

if nargin < 3, x = 0; else x = x-1; end
if nargin < 4, y = 0; else y = y-1; end
if nargin < 5, h = r.getSizeX(); end
if nargin < 6, w = r.getSizeY(); end
   
plane = r.openBytes(iPlane - 1, x,y,h,w);
    
% convert byte array to MATLAB image
if sgn
    % can get the data directly to a matrix
    I = loci.common.DataTools.makeDataArray2D(plane, ...
        bpp, fp, little, ip.Results.height);
else
    % get the data as a vector, either because makeDataArray2D
    % is not available, or we need a vector for typecast
    I = loci.common.DataTools.makeDataArray(plane, ...
        bpp, fp, little);
end

% Java does not have explicitly unsigned data types;
% hence, we must inform MATLAB when the data is unsigned
if ~sgn
    % NB: arr will always be a vector here
    switch class(I)
        case 'int8'
            I = typecast(I, 'uint8');
        case 'int16'
            I = typecast(I, 'uint16');
        case 'int32'
            I = typecast(I, 'uint32');
        case 'int64'
            I = typecast(I, 'uint64');
    end
end

if isvector(I)
    % convert results from vector to matrix
    shape = [h, w];
    I = reshape(I, shape)';
end

end



