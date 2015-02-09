function [data, metatdata, ireader] = imreadBF(name, varargin) 
%
% [data, metatdata] = imreadBF(name, param) 
%
% description:
%    import bio-format file using the loci tools
%
% input:
%    name   (optional)    filename, if not gived a file borwser to select file is opened
%    param  (optional)    struct with optional entries
%           .S / .s       ids of series to import (if non scalar cell array is returned) ([] = all)
%           .T / .t       ids of time frames to import [tid1, tid2, ...] or {tid1,tid1, ...} ([] = all)
%           .C / .c       ids of channels to import [cid1, cid2, ...] or {cid1, cid2,...} ([] = all)
%           .Z / .z       pixel ids in z direction [zid1, zid2, ...] or {zid1,zid1, ...} ([] = all)
%           .Y / .y       pixel ids in y direction [yid1, yid2, ...] or {yid1,yid1, ...} ([] = all)
%           .X / .x       pixel ids in x direction [xid1, xid2, ...] or {xid1,xid1, ...} ([] = all)
%           .metadata     read meta data (true)
%           .maxmemory    maximal estimated size of image in memory in GB before reading is canceled (5)
%           .squeeze      squeeze the returned data (true)
%           .gui          open gui to select series if there are multiple (full image series is imported by default)
%           .verbose      print some info when reading images (true)
%           .close        colse the passed reader
%
% output:
%    data     image data
%    metadata (optional) image meta data 
%    ireader  (optional) image reader from loci tools for further image access
%
% note: 
%    the loci_tools.jar must be present in the javaclasspath
%    try bfinitialize or ijintialize to set this path
%
% See also: bfinitialize, imread, ijimage2mat

param = parseParameter(varargin);

gui = getParameter(param, 'gui', false);

if gui
   data = importlocigui(name);
else 
   ireader = imreadBFReader(name);
   [data, metatdata] = importlociseries(ireader, param);   

   if nargout < 3 && getParameter(param, 'close', ~isa(name, 'loci.formats.ChannelSeparator'));
      ireader.close();
   end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main data import function 

function [data, metadata] = importlociseries(r, param)

   % determine sizes
   N = r.getSeriesCount();
   T = r.getSizeT();
   C = r.getSizeC();
   Z = r.getSizeZ();
   Y = r.getSizeY();
   X = r.getSizeX();

   % determine images to read
   s = imfrmtParseRangeToIndex(N, 'S', param);

   % read multiple series into cells
   if length(s) > 1
      par = param;
      if isfield(par, 'S')
         par = rmfield(par, 'S');
      end
      if isfield(par, 's')
         par = rmfield(par, 's');
      end
      
      ns = length(s);
      data = cell(ns,1);
      metadata = cell(ns,1);
      for i = 1:length(s)
         par.S = s(i);
         [d, m] = importlociseries(r, par);
         data{i} = d; metadata{i} = m;
      end
      return
   end

   %%% read a single series 
   if s < 1 || s > N
      error('imreadBF: series index %g out of range %g : %g!', s, 1, N);
   end

   r.setSeries(s - 1);
   n = r.getImageCount();
 
   % dermine ids to read
   tids = imfrmtParseRangeToIndex(T, 'T', param);
   cids = imfrmtParseRangeToIndex(C, 'C', param);
   zids = imfrmtParseRangeToIndex(Z, 'Z', param);
   yids = imfrmtParseRangeToIndex(Y, 'Y', param);
   yids = Y - yids + 1;  % y direction in bf from top to bottom thus inverse for Y
   xids = imfrmtParseRangeToIndex(X, 'X', param);

   % determine start and end ids for x and y
   xs = min(xids); 
   ys = min(yids); 
   
   w = max(xids)-min(xids)+1;
   h = max(yids)-min(yids)+1;

   xids = xids - xs + 1;
   yids = yids - ys + 1;

   % final sizes after ids selection
   T = length(tids);
   C = length(cids);
   Z = length(zids);
   Y = length(yids);
   X = length(xids);

   % check memory constraints
   ps = loci.formats.FormatTools.getPlaneSize(r, h, w);
   mm = getParameter(param, 'maxmemory', 5);
   if ps/(1024)^3 >= mm
       error('imreadBF: image larger than %gGB! increase maxmemory parameter!', mm);
   end

   % print some info
   verbose = getParameter(param, 'verbose', true);
   if verbose 
      fprintf('imreadBF: series %d, size: %g x %g x %g x %g x %g\n', s, X, Y, Z, C, T);
   end
   printprogress  = verbose && n > 5;
  
   % allocate memory
   data = zeros([X, Y, Z, C, T]);

   % pixel types, bytes per pixel, etc
   pixelType     = r.getPixelType();
   bytesPerPixel = loci.formats.FormatTools.getBytesPerPixel(pixelType);
   
   isSigned        = loci.formats.FormatTools.isSigned(pixelType);
   isLittleEndian  = r.isLittleEndian();
   isFloatingPoint = loci.formats.FormatTools.isFloatingPoint(pixelType);
   

   % loop over all images
   for i = 1:n 
      
      % print progress
      if printprogress
         if mod(i+1, 72) == 1
            fprintf('\n    ');
         end
         fprintf('.');
      end

      % determine coordinates 
      zct = r.getZCTCoords(i-1) + 1;
      z = find(zids == zct(1), 1, 'first');
      c = find(cids == zct(2), 1, 'first');
      t = find(tids == zct(3), 1, 'first');
      
      % color maps data
      %if bytesPerPixel == 1
      %  colorMaps{i} = r.get8BitLookupTable()';
      %else
      %  colorMaps{i} = r.get16BitLookupTable()';
      %end

      % read if required
      if ~isempty(z) && ~isempty(c) && ~isempty(t)
         
         % load image plane
         %fprintf('imreadBF: reading plane: %g %g %g\n', z, c, t);
         plane = r.openBytes(i - 1, xs - 1, ys - 1, w, h);
         img = loci.common.DataTools.makeDataArray(plane, bytesPerPixel, isFloatingPoint, isLittleEndian);

         % cast to data type
         if ~isSigned
            switch class(img)
               case 'int8'
                  img = typecast(img, 'uint8');
               case 'int16'
                  img = typecast(img, 'uint16');
               case 'int32'
                  img = typecast(img, 'uint32');
               case 'int64'
                  img = typecast(img, 'uint64');
            end
         end

         % reshape to pq format
         if isvector(img)
            img = reshape(img, [w, h]);
         end
         
         %assign data
         data(:,:,z,c,t) = img(xids, yids);
      end
   end

   if printprogress
      fprintf('\n');
   end

   % squeeze
   if getParameter(param, 'squeeze', true)
      data = squeeze(data);
   end

   % get meta data if required
   if getParameter(param, 'metadata', true)
      metadata = r.getSeriesMetadata();
   else
      metadata = [];
   end
   
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GUI interface to open image via seris selection
function data = importlocigui(name)
frames
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







    
 