function [data, metatdata] = imread_bf(name, varargin) 
%
% [data, metatdata] = imread_bf(name, param) 
%
% description:
%    import bio-format file by opening the loci importer gui
%
% input:
%    name   (optional)   filename
%    param  (optional)   struct with optional entries (if not present gui is opened, [] imports everything)
%           .series      ids of series to import (if array data is cell array) ([] = all)
%           .time /.t    ids of time frames to import [tid1, tid2, ...] or {tid1,tid1, ...} ([] = all)
%           .channel /.c ids of channels to import [cid1, cid2, ...] or {cid1, cid2,...} ([] = all)
%           .z /.l       pixel ids in z / l direction [lid1, lid2, ...] or {lid1,lid1, ...} ([] = all)
%           .x /.p       pixel ids in x / p direction [pid1, pid2, ...] or {pid1,pid1, ...} ([] = all)
%           .y /.q       pixel ids in y / q direction [qid1, qid2, ...] or {qid1,qid1, ...}{ymin, ymax} ([] = all)
%           .metadata    read meta data (true)
%           .gui         open gui to select series if there are multiple (full image series is imported by default)
%           .maxmem      maximal estimated size of image in memory in GB (5)
%           .squeeze     squeeze the returned data (true)
%
% output:
%    data     image data
%    metadata (optional) image meta data 
%
% note: the loci_tools.jar must be present in the javaclasspath
%       try bfinitialize or ijintialize to set this path
%
% See also: bfinitialize, imread, ijimage2mat


% Prompt for a file if no input
if nargin == 0
  [file, path] = uigetfile(bffileextensions, 'Choose a file to open');
  name = [path file];
end
if ~exist(name, 'file')
   error('imread_bf: file does not exists: %s', name)
end

param = parseParameter(varargin{:});

% check for multiple series
% if nargin < 2
%    r = loci.formats.ChannelFiller();
%    r = loci.formats.ChannelSeparator(r);
%    r.setId(name);
%    numSeries = r.getSeriesCount();
%  
%    if numSeries > 1
%       param.gui = true;
%    else
%       param.gui = false;
%    end
% end

gui = getParameter(param, 'gui', false);

if gui
   data = importlocigui(name);
else
   [data, metatdata] = importlociseries(name, param);
end

end



function [data, metadata] = importlociseries(name, param)

% initialize logging
% loci.common.DebugTools.enableLogging('INFO');

% Get the channel filler
r = loci.formats.ChannelFiller();
r = loci.formats.ChannelSeparator(r);
%if stitchFiles
%    r = loci.formats.FileStitcher(r);
%end
r.setMetadataStore(loci.formats.MetadataTools.createOMEXMLMetadata());
r.setId(name);

numSeries = r.getSeriesCount();

s = getParameter(param, 'series', []);

if isempty(s)
   s = 1:numSeries;
end

if ~iscell(s) && numel(s) > 1
   s = num2cell(s(:));
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
   error('imread_bf: series index %g out of range %g : %g!', s, 1, numSeries);
end

%%% read a single series 

r.setSeries(s - 1);
numImages = r.getImageCount();

% dimensions
sizeT = r.getSizeT();
sizeC = r.getSizeC();
sizeZ = r.getSizeZ();
sizeY = r.getSizeY();
sizeX = r.getSizeX();
%var2char({sizeX, sizeY, sizeZ, sizeC, sizeT})

tids = getIds(param, 'time'   , 't', sizeT);
cids = getIds(param, 'channel', 'c', sizeC);

% z is reverse l !
zids = getIds(param, 'z',       'l', sizeZ);
if isentry(param, 'z')
   zids = sizeZ - zids + 1;
end

% y is reverse q !
yids = getIds(param, 'y',       'q', sizeY);
if ~isentry(param, 'y')
   yids = sizeY - yids + 1;
end
xids = getIds(param, 'x',       'p', sizeX);


%determine start and end ids for x and y
xs = min(xids); w = max(xids)-min(xids)+1;
ys = min(yids); h = max(yids)-min(yids)+1;

xids = xids - xs + 1;
yids = yids - ys + 1;

sizeT = length(tids);
sizeC = length(cids);
sizeZ = length(zids);
sizeY = length(yids);
sizeX = length(xids);

% enogh memory -> get plane size
planeSize = loci.formats.FormatTools.getPlaneSize(r, h, w);
%planeSize = loci.formats.FormatTools.getPlaneSize(r);

if planeSize/(1024)^3 >= getParameter(param, 'maxmem', 5)
    error('imread_bf: image larger than %gGB! try opening images individually !', getParameter(param, 'maxmem', 5));
end

fprintf('imread_bf: size of image: %g x %g x %g x %g x %g\n', sizeX, sizeY, sizeZ, sizeC, sizeT);
data = zeros([sizeX, sizeY, sizeZ, sizeC, sizeT]);


printprogress  = numImages > 5;
if printprogress
   fprintf('Reading series #%d', s);
end
for i = 1:numImages
   if printprogress
      if mod(i+1, 72) == 1
         fprintf('\n    ');
      end
      fprintf('.');
   end
      
   zct = r.getZCTCoords(i-1) + 1;
   z = find(zids == zct(1), 1, 'first');
   c = find(cids == zct(2), 1, 'first');
   t = find(tids == zct(3), 1, 'first');
  
   
   if ~isempty(z) && ~isempty(c) && ~isempty(t)
      %fprintf('%g %g %g\n', z, c, t);
      %size( bfgetplane(r, i, x, y, h, w))
      bfdata = bfgetplane(r, i, xs, ys, w, h);   
      
      %size(bfdata)
      
      data(:,:,z,c,t) = bfdata(xids, yids);
   end 
end

if printprogress
   fprintf('\n');
end

if getParameter(param, 'squeeze', true)
   data = squeeze(data);
end

if getParameter(param, 'metadata', true)
   metadata = r.getSeriesMetadata();
else
   metadata = [];
end

r.close();

end


function I = bfgetplane(r, iPlane, x, y, w, h)
   % Get pixel type
   pixelType = r.getPixelType();
   bpp = loci.formats.FormatTools.getBytesPerPixel(pixelType);
   fp = loci.formats.FormatTools.isFloatingPoint(pixelType);
   sgn = loci.formats.FormatTools.isSigned(pixelType);
   little = r.isLittleEndian();

%    if nargin < 3, x = 0; else x = x-1; end
%    if nargin < 4, y = 0; else y = y-1; end
%    if nargin < 5, h = r.getSizeX(); end
%    if nargin < 6, w = r.getSizeY(); end

   plane = r.openBytes(iPlane - 1, x - 1, y - 1, w, h);

   % convert byte array to MATLAB image
   if sgn
      % can get the data as matrix
      I = loci.common.DataTools.makeDataArray2D(plane, bpp, fp, little, ip.Results.height);
   else
      % get the data as a vector for typecast
      I = loci.common.DataTools.makeDataArray(plane, bpp, fp, little);
   end

   % Java does not have explicitly unsigned data types;
   % hence, we must inform MATLAB when the data is unsigned
   if ~sgn
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

   % convert results from vector to matrix in pq format
   if isvector(I)
      I = reshape(I, [w, h]);
   end
   
end


function ids = getIds(param, name, nameshort, maxid)
      
   if isentry(param, name)
      ids = getParameter(param, name, []);
   else
      ids = getParameter(param, nameshort, []);
   end
   

   if isempty(ids)
      ids = 1:maxid;
   elseif iscell(ids)
      %if length(ids) ~= 2
      %   error(['imread_bf: ', name, ' indices not a cell of {minid,maxid} but %s!'], var2char(ids))
      %end
      
      %if isempty(ids{2}) || ischar(ids{2}) || ids{2} > maxid  % char = 'end'
      %   ids{2} = maxid;
      %end
      %if isempty(ids{1}) || ischar(ids{1}) || ids{1} < 1
      %   ids{1} = 1;
      %end
      %ids = ids{1} : ids{2};
      
      ids = cell2mat(ids);
   end

   ids = ids(:)';
   
   if max(ids) > maxid || min(ids) < 1
      error(['imread_bf: ', name, ' indices out of range !'])
   end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GUI interface to open image via seris selection
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
         data{i} = ijimage2mat(imps(i));
      end
   end

end






%%% useful code


%pixelType = r.getPixelType();
%bpp = loci.formats.FormatTools.getBytesPerPixel(pixelType);
%bppMax = power(2, bpp * 8);
%imageList = cell(numImages, 2);
%colorMaps = cell(numImages);


%arr = bfGetPlane(r, i, varargin{:});

%   % retrieve color map data
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
   
   
%     % build an informative title for our figure
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


% save images and metadata
%result{s, 1} = imageList;
%result{s, 2} = r.getSeriesMetadata();
%result{s, 3} = colorMaps;
%result{s, 4} = r.getMetadataStore();
