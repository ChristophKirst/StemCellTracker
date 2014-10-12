classdef ImageInfo < matlab.mixin.Copyable
   
   % class that manages information about a single image in a file
   % for multiple series in a file choose ImageSourceList
   
   % notes:
   %    the info class is based on the info needed to read images via imreadBF / loci tools
   %    genreal organization:
   %                the ouput is represented in the data and cell size/format information
   %                raw input formats are stored in the raw data dn cell size / formats
   %                and possible reshaping of the raw cell data is done via the raw cell reshape options
   %                data flow: read image gives raw data when cell ranges are specified to be a singelton
   %                           the cell organizaton for reading is in the raw cell format
   %                           the cell data is then reshaped and reformatted into the output data and cells
   %                           with possible conversions between cell and data dimensions
   
   properties
      % output data/cell sizes and formats when returned by the data and cell commands 
      idatasize      = [];           % size of the image 
      idataformat    = '';           % format of the image 

      idataclass = '';               % class of the image 
      
      icellsize   = 1;               % size of the cell structure
      icellformat = '';              % format of the cell structure
      
      irange = struct();             % specify a subset of coordinate ranges
 
      % input  data/cell sizes and formats when load from disk / stream etc
      irawdatasize = [];             % size of the raw data, [] = idatasize
      irawdataformat = '';           % format of raw data obtained by reading routines
      
      irawdataclass = '';            % class of the image 
      
      irawcellsize = [];             % size of the raw cell structure
      irawcellformat = '';           % format of the raw cells structure
      
      % reshaping
      ireshapefrom = {};             % formats dims in raw for reshaping 
      ireshapeto   = {};             % formats before and after reshaping each raw data array
      ireshapesize = {};             % sizes for the reshaping
           
      % meta
      ichannelname = [];             % channel name
      icolor = {};                   % colors for the color channel C/c
      icoloralpha = [];              % (optional) alpha values for the color channel C/c (=1 if no specified) 
  
      % intensity scale
      imaxdataintensity = [];        % minimal intensity of data
      imindataintensity = [];        % maximal intensity of data
      
      imaxrawintensity = [];         % minimal intensity of raw data
      iminrawintensity = [];         % maximal intensity of raw data
      
      % spatial scales
      iscale = containers.Map;       % (optional) spatial scale of image in distance/time per pixel (coords: X, Y, Z, T)
      iunit  = containers.Map;       % (optional) spatial unit of image in distance/time per pixel (coords: X, Y, Z, T)
      
      iname  = [];                   % (optional) name of the image  
              
      imetadata = [];                % (optional) meta data
   end
   
   methods

      function obj = ImageInfo(varargin) % constructor
         %
         % ImageInfo()
         % ImageInfo(data)
         % ImageInfo(...,fieldname, fieldvalue,...)
         %
         
         if nargin == 0
            return
         elseif nargin == 1
            if isa(varargin{1}, 'ImageInfo') %% copy constructor
               obj = copy(varargin{1});
            elseif isnumeric(varargin{1})
               obj.fromData(varargin{1});
            elseif ischar(varargin{1})
               obj.fromFile(varargin{1});
            else
               error('%s: invalid constructor input, expects char at position %g',class(obj), 1);
            end
         else
            obj = obj.fromParameter(varargin);
         end
      end
      
      function obj = fromParameter(obj, varargin)
         disp info
         obj = classFromParameter(obj, 'i', varargin);
         obj.initializeRawFormatsAndSizes();
      end
      
      function obj = fromData(obj, data)
         obj = imfrmtInfoFromData(data, obj);
         obj.initializeRawCellDataFormatsAndSizesFromData();
      end

      function obj = fromFile(obj, filename)
         info = imreadBFInfo(filename);
         for n = properties(obj)'
            obj.(n{1}) = info.(n{1});
         end
         obj.initializeRawCellDataFormatsAndSizesFromData();
      end
      
      function obj = fromImageInfo(obj, info)
         for n = properties(ImageInfo)'
            obj.(n{1}) = info.(n{1});
         end
      end
      
      function obj = fromImfInfo(obj, info)
         info = imfinfo2info(info);
         for n = properties(obj)'
            obj.(n{1}) = info.(n{1});
         end
         obj.initializeRawCellDataFormatsAndSizesFromData();
      end

      function obj = initializeRawCellDataFormatsAndSizesFromData(obj)
         if isempty(obj.irawdataformat)
            obj.irawdataformat = obj.idataformat;
            obj.irawdatasize   = obj.idatasize;
         end
         
         if isempty(obj.irawcellformat)
            obj.irawcellformat = obj.icellformat;
            obj.irawcellsize   = obj.icellsize;
         end
      end
      
      function obj = initializeCellDataFormatsAndSizesFromRaw(obj)
         if isempty(obj.idataformat)
            obj.idataformat = obj.irawdataformat;
            obj.idatasize   = obj.irawdatasize;
         end
         
         if isempty(obj.icellformat)
            obj.icellformat = obj.irawcellformat;
            obj.icellsize   = obj.irawcellsize;
         end
      end
      
 
      function obj = renameFormat(obj, oldlab, newlab)
         obj.renameDataFormat(oldlab, newlab);
         obj.renameRawDataFormat(oldlab, newlab);
        
         obj.renameCellFormat(oldlab, newlab);
         obj.renameRawCellFormat(oldlab, newlab);
         
         obj.irange = renamestruct(obj.irange, num2cell(oldlab), num2cell(newlab));
      end
      
      function obj = initializeDataClassFromRaw(obj)
         obj.idataclass = obj.irawdataclass;
      end
      
            
      function obj = initializeDataAndCellSizeFromRaw(obj)
         %
         % obj = initializeDataAndCellSizeFromRaw(obj)
         %
         % decription:
         %   sets data and cell size using raw data and cell size and format info 
         
         % check if dimensions are lost and if so set iranges to first 
         outfrmt= imfrmtReshapeFormat(obj.rawCellDataFormat, obj.reshapeFrom, obj.reshapeTo);
         cdfrmt = obj.cellDataFormat;
         
         extrafrmt = imfrmtExtraFormats(cdfrmt, outfrmt);

         if ~isempty(extrafrmt)   
            % size of extra dims
            extrasize = imfrmtReshapeCellDataSize( obj.rawDataSize, obj.rawCellSize, obj.rawDataFormat, obj.rawCellFormat, ...
                                                   extrafrmt, '', obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize);
                                              
            if any(extrasize > 1)
               warning('%s: input dimensions %s of size %s lost and set to singelton!', class(obj), extrafrmt, var2char(extrasize))
            
               
               rgs = [num2cell(extrafrmt); num2cell(ones(1,length(extrafrmt)))];
               rgs = imfrmtParseRangeLast(struct(rgs{:}), obj.irange);
               rgs = imfrmtParseRange(extrasize, extrafrmt, rgs);
               rgsn = fieldnames(rgs);
               for i = 1:length(rgsn)
                  v= rgs.(rgsn{i});
                  rgs.(rgsn{i}) = v(1);
               end
               obj.irange = rgs;
            end
         end
 
         % reshaped sizes
         [obj.idatasize, obj.icellsize] = imfrmtReshapeCellDataSize( obj.rawDataSize, obj.rawCellSize, obj.rawDataFormat, obj.rawCellFormat, ...
                                                                     obj.dataFormat, obj.cellFormat, obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize);       
              
         % correct for ranges                                                             
         if ~isempty(obj.irange)
            obj.idatasize = imfrmtRangeSize(obj.dataSize, obj.dataFormat, obj.irange);
            obj.icellsize = imfrmtRangeSize(obj.cellSize, obj.cellFormat, obj.irange);
         end
      end
      
      function obj = initializeDataAndCellFormatFromRaw(obj)
         obj.idataformat = obj.irawdataformat;
         obj.icellformat = obj.irawcellformat;
      end
      
      function obj = initializeDataAndCellFromRaw(obj)
         obj.initializeDataAndCellFormatFromRaw();
         obj.initializeDataAndCellSizeFromRaw();
         obj.initializeDataClassFromRaw();
      end

      
      function obj = initializeReshape(obj)
         if ~iscell(obj.ireshapefrom)
            obj.ireshapefrom = {obj.ireshapefrom};
         end
         
         if ~iscell(obj.ireshapeto)
            obj.ireshapeto = {obj.ireshapeto};
         end
         
         if ~iscell(obj.ireshapesize)
            obj.ireshapesize = {obj.ireshapesize};
         end

         if length(obj.ireshapesize) ~= length(obj.ireshapefrom) || length(obj.ireshapesize) ~= length(obj.ireshapeto)
            error('ImageInfo: initializeReshape: inconsistent reshape formats or sizes');
         end
         
         [obj.idataformat, obj.icellformat] = imfrmtReshapeCellDataFormat(obj.irawdataformat, obj.irawcellformat, obj.ireshapefrom, obj.ireshapeto);

         obj.initializeDataAndCellSizeFromRaw;
      end

      function obj = setReshape(obj, reshapeFrom, reshapeTo, reshapeSize)
         obj.ireshapefrom = reshapeFrom;
         obj.ireshapeto   = reshapeTo;
         obj.ireshapesize = reshapeSize;   
         obj.initializeReshape();
      end

      function obj = addReshape(obj, reshapeFrom, reshapeTo, reshapeSize)
         obj.ireshapefrom{end+1} = reshapeFrom;
         obj.ireshapeto{end+1}   = reshapeTo;
         obj.ireshapeSize{end+1} = reshapeSize;
         obj.initializeReshape();
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% data
      
      function d = dataDims(obj)
         d = length(obj.idataformat);
      end 

      function s = dataSize(obj, varargin)
         if nargin == 1
            s = obj.idatasize;
         else
            s = imfrmtReformatSize(obj.idatasize, obj.idataformat, varargin{1});
         end
      end
      
      function ds = dataSizeXYZCT(obj)
         ds = imfrmtReformatSize(obj.idatasize, obj.idataformat, 'XYZCT');
      end

      function ds = dataSizeX(obj)
         ds = imfrmtReformatSize(obj.idatasize, obj.idataformat, 'X');
      end
      function ds = dataSizeY(obj)
         ds = imfrmtReformatSize(obj.idatasize, obj.idataformat, 'Y');
      end
      function ds = dataSizeZ(obj)
         ds = imfrmtReformatSize(obj.idatasize, obj.idataformat, 'Z');
      end
      function ds = dataSizeC(obj)
         ds = imfrmtReformatSize(obj.idatasize, obj.idataformat, 'C');
      end
      function ds = dataSizeT(obj)
         ds = imfrmtReformatSize(obj.idatasize, obj.idataformat, 'T');
      end

      function df = dataFormat(obj)
         df = obj.idataformat;
      end
      
      function dc = dataClass(obj)
         dc = obj.idataclass;
      end
      
      
      function obj = setDataSize(obj, newsize)
         obj.idatasize = newsize;
      end
      
      function obj = setDataFormat(obj, newfrmt)
         obj.idataformat = newfrmt;  
      end
      
      function obj = setDataFormatAndSize(obj, newfrmt, newsize)
         obj.idatasize = newsize;
         obj.idataformat = newfrmt;
      end

      function obj = reformatDataFormat(obj, newfrmt)
         obj.reformatColor(obj.idataformat, newfrmt);
         obj.idatasize   = imfrmtReformatSize(obj, obj.idatasize, obj.idatafrmt, newfrmt);
         obj.idataformat = newfrmt;
      end

      function obj = renameDataFormat(obj, oldlab, newlab) 
          %
          % obj = renameDataFormat(obj, oldlab, newlab) 
          %
          % description:
          %    exchanges names of the dimensions oldlab(i) -> newlab(i)
          %    
          % input:
          %    oldlab     the format dimensions to change
          %    newlab     the new format dimensions 
           
          % check for conflict
          if any(ismember(setdiff(obj.idataformat, oldlab), newlab)) || length(oldlab) ~= length(newlab)
             error('%s: renameDataFormat: inconsistent reformatting of format %s from %s to %s', class(obj), obj.idataformat, oldlab, newlab)
          end    
          % rename
          ids = ismember(obj.idataformat, oldlab);
          obj.idataformat(ids) = newlab;
          
          obj.irange = renamestruct(obj.irange, num2cell(oldlab), num2cell(newlab));
      end
       
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% cell

      function d = cellDim(obj)
         d = length(obj.icellformat);
      end

      function f = cellFormat(obj)
         f = obj.icellformat;
      end
      
      function s = cellSize(obj)
         s = obj.icellsize;
      end
      
      function obj = setCellSize(obj, newsize)
         obj.isize = newsize;
      end
         
      function obj = setCellFormat(obj, newfrmt)
         obj.icellformat = newfrmt;
      end

      function obj = reformatCellFormat(obj, newfrmt)
         obj.icellsize  = imfrmtReformatSize(obj, obj.cellSize, obj.cellFormat, newfrmt);
         obj.icellformat = newfrmt;
      end 
      
      function obj = renameCellFormat(obj, oldlab, newlab)            
          % check for conflict
          if any(ismember(setdiff(obj.icellformat, oldlab), newlab)) || length(oldlab) ~= length(newlab)
             error('%s: renameRawFormat: inconsistent reformatting of format %s from %s to %s', class(obj), obj.icellformat, oldlab, newlab)
          end    
          % rename
          ids = ismember(obj.icellformat, oldlab);
          obj.icellformat(ids) = newlab;
          
          obj.irange = renamestruct(obj.irange, num2cell(oldlab), num2cell(newlab));
      end

      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% raw data
        
 
      function d = rawDataDims(obj)
         d = length(obj.rawDataFormat);
      end 

      function rs = rawDataSize(obj, varargin)
         if nargin == 1
            rs = obj.irawdatasize;
         else
            rs = imfrmtReformatSize(obj.irawdatasize, obj.irawdataformat, varargin{1});
         end
      end
      
      function rf = rawDataFormat(obj)
         rf = obj.irawdataformat;
      end
      
      function rc = rawDataClass(obj)
         rc = obj.irawdataclass;
      end
 
      function obj = setRawDataSize(obj, newsize)
         obj.irawdatasize = newsize;
      end
      
      function obj = setRawDataFormat(obj, newfrmt)
         obj.irawdataformat = newfrmt;
      end
      
      function obj = reformatRawDataFormat(obj, newfrmt)
         obj.irawdatasize   = imfrmtReformatSize(obj.rawDataSize, obj.rawDataFormat, newfrmt);
         obj.irawdataformat = newfrmt;
      end
      
      function obj = renameRawDataFormat(obj, oldlab, newlab)            
          % check for conflict
          if any(ismember(setdiff(obj.irawdataformat, oldlab), newlab)) || length(oldlab) ~= length(newlab)
             error('%s: renameRawFormat: inconsistent reformatting of format %s from %s to %s', class(obj), obj.irawdataformat, oldlab, newlab)
          end    
          % rename
          ids = ismember(obj.irawdataformat, oldlab);
          obj.irawdataformat(ids) = newlab;
      end

      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% raw cell

      function d = rawCellDims(obj)
         d = length(obj.rawCellFormat);
      end 
      
      function rs = rawCellSize(obj, varargin)
         if nargin == 1
            rs = obj.irawcellsize;
         else
            rs = imfrmtReformatSize(obj.irawcellsize, obj.irawcellforamt, varargin{1});
         end
      end
      
      function rf = rawCellFormat(obj)
         rf = obj.irawcellformat;
      end
      
      function obj = setRawCellSize(obj, newsize)
         obj.irawcellsize = newsize;
      end
      
      function obj = setRawCellFormat(obj, newfrmt)
         obj.irawcellformat = newfrmt;
      end
      
      function obj = reformatRawCellFormat(obj, newfrmt)
         obj.irawcellsize  = imfrmtReformatSize(obj, obj.rawCellSize, obj.rawCellFormat, newfrmt);
         obj.irawcellformat = newfrmt;
      end
      
      function obj = renameRawCellFormat(obj, oldlab, newlab)            
          % check for conflict
          if any(ismember(setdiff(obj.irawcellformat, oldlab), newlab)) || length(oldlab) ~= length(newlab)
             error('%s: renameRawFormat: inconsistent reformatting of format %s from %s to %s', class(obj), obj.irawcellformat, oldlab, newlab)
          end    
          % rename
          ids = ismember(obj.irawcellformat, oldlab);
          obj.irawcellformat(ids) = newlab;
      end
      
      
        
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% celldata / full format 
      
      function f = cellDataFormat(obj)
         f = [obj.dataFormat, obj.cellFormat];
      end

      function f = cellDataSize(obj)
         f = [obj.dataSize, obj.dataFormat];
      end
          
      function f =rawCellDataFormat(obj)
         f = [obj.rawDataFormat, obj.rawCellFormat];
      end

      function f = rawCellDataSize(obj)
         f = [obj.rawDataSize, obj.rawDataFormat];
      end
      
      
      function obj = setCellDataFormat(obj, newdatafrmt, newcellfrmt)
         obj.idataformat = newdatafrmt;
         obj.icellformat = newcellfrmt;
         obj.initializeDataAndCellSizeFromRaw();
      end

      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% reshaping 
          
      function r = reshapeFrom(obj)
         r = obj.ireshapefrom;
      end
      
      function r = reshapeTo(obj)
         r = obj.ireshapeto;
      end
      
      function r = reshapeSize(obj)
         r = obj.ireshapesize;
      end
        
      function df = dataFormatFromRaw(obj, varargin)
         % uses reshaping and raw format info to get version of the data frmt
         df = imfrmtReshapeFormat(obj.rawDataFormat, obj.reshapeFrom. obj.reshapeTo); 
      end
      
      function cf = cellFormatFromRaw(obj, varargin)
         % uses reshaping and raw format info to get version of the cell frmt
         cf = imfrmtReshapeFormat(obj.rawCellFormat, obj.reshapeFrom. obj.reshapeTo);  
      end

      
      function d = reshapeData(obj, d, varargin)
         % uses reshaping and format info to reshape the raw data d to the return data
         d = imfrmtReshape(d, obj.rawDataFormat, obj.dataFormat, obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize);  
      end
      
      function d = reshapeCell(obj, d, varargin)
         % uses reshaping and format info to reshape the raw data d to the return data
         d = imfrmtReshape(d, obj.rawCellFormat, obj.cellFormat, obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize);  
      end

      function d = reshapeCellData(obj, d, varargin)
         % uses reshaping and format info to reshape the raw data d to the return data
         d = imfrmtReshapeCellData(d,obj.rawDataFormat,  obj.rawCellFormat, obj.dataFormat, obj.cellFormat, obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize);  
      end

      function [rawRange, rawReshape] = reshapeRangeToRawRange(obj, varargin)
         % uses reshaping and format info to reshape the raw data d to the return data
         rawRange = imfrmtParseRangeLast(obj.irange, varargin);
         [rawRange, rawReshape] = imfrmtReshapeInverseCellDataRange(obj.dataSize, obj.dataFormat, obj.rawDataFormat, obj.cellSize, obj.cellFormat, obj.rawCellFormat, ...
                                                    obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize, rawRange);
      end
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% ranges
          
      function r = range(obj)
         r = obj.irange;
      end
      
      function rn = rangeNames(obj)
         rn = fieldnames(obj.irange);
      end
      
      function obj = setRange(obj, varargin)
         obj.irange = imfrmtParseRange(obj.cellDataSize, obj.cellDataFormat, varargin);
         obj.initializeDataAndCellSizeFromRaw;
      end
      
      function obj = addRange(obj, name, range)
         obj.irange.(name) = range;
         obj.initializeDataAndCellSizeFromRaw;
      end
      
      function obj = resetRange(obj)
         obj.irange = struct();
         obj.initializeDataAndCellSizeFromRaw;
      end
      
      
      
      

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% assignments
      
      function asgn = dataAssignment(obj, varargin)
         asgn = imfrmtAssignment(obj.dataSize, obj.dataFormat, varargin); 
      end
      
      function asgn = cellAssignment(obj, varargin)
         asgn = imfrmtAssignment(obj.cellSize, obj.cellFormat, varargin);
      end
      
      function asgn = rawDataAssignment(obj, varargin)
         asgn = imfrmtAssignment(obj.rawDataSize, obj.rawDataFormat, varargin); 
      end
      
      function asgn = rawCellAssignment(obj, varargin)
         asgn = imfrmtAssignment(obj.rawCellSize, obj.rawCellFormat, varargin);
      end
      
      
      function i = dataFormatPosition(obj, frmt)
         %
         % i = dataFormatPosition(obj, frmt)
         % 
         % description: 
         %    returns position of the format names in frmt in obj.dataFormat
         
         [~, i] = ismember(frmt, obj.dataFormat);
      end
      
      function i = cellFormatPosition(obj, frmt)
         [~, i] = ismember(frmt, obj.cellFormat);
      end
 
      function i = rawDataFormatPosition(obj, frmt)
         [~, i] = ismember(frmt, obj.rawDataFormat);
      end
      
      function i = rawCellFormatPosition(obj, frmt)
         [~, i] = ismember(frmt, obj.rawCellFormat);
      end
  
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% indices
      
            
      function n = nDataIndices(obj)
         n = prod(obj.dataSize);
      end
      
      function n = nCellIndices(obj)
         n = prod(obj.cellSize);
      end
      
      function n = nCells(obj)
         n = prod(obj.cellSize);
      end
      
      function n = nRawDataIndices(obj)
         n = prod(obj.rawDataSize);
      end
      
      function n = nRawCellIndices(obj)
         n = prod(obj.rawCellSize);
      end

      function i = dataIndex(obj, varargin)
         if nargin >= 2 && isnumeric(varargin{1}) && isscalar(varargin{1})
            if nargin > 2
               i = sub2ind(obj.dataSize, varargin{:});
            else
               i = varargin{1};  % index list or array
            end
         else
            i = imfrmtRangeToIndex(obj.dataSize, obj.dataFormat, parseParameter(varargin));
         end
      end
      
      function i = cellIndex(obj, varargin)
         if nargin >= 2 && isnumeric(varargin{1})
            if nargin > 2
               i = sub2ind(obj.dataSize, varargin{:});
            else
               i = varargin{1};
            end
         else
            i = imfrmtRangeToIndex(obj.cellSize, obj.cellFormat, parseParameter(varargin));
         end
      end
      
      function i = rawDataIndex(obj, varargin)
         if nargin >= 2 && isnumeric(varargin{1})
            if nargin > 2
               i = sub2ind(obj.dataSize, varargin{:});
            else
               i = varargin{1};
            end
         else
            i = imfrmtRangeToIndex(obj.rawDataSize, obj.rawDataFormat, parseParameter(varargin));
         end
      end
      
      function i = rawCellIndex(obj, varargin)
         if nargin >= 2 && isnumeric(varargin{1})
            if nargin > 2
               i = sub2ind(obj.dataSize, varargin{:});
            else
               i = varargin{1};
            end
         else
            i = imfrmtRangeToIndex(obj.rawCellSize, obj.rawCellFormat, parseParameter(varargin));
         end
      end
      
      
      function coords = dataIndexToCoordinate(obj, id, varargin)
         coords = imfrmtIndexToCoordinate(obj.dataSize, obj.dataFormat, id);
      end
      
      function coords = cellIndexToCoordinate(obj, id, varargin)
         coords = imfrmtIndexToCoordinate(obj.cellSize, obj.cellFormat, id);
      end
      
      
      function coords = rawDataIndexToCoordinate(obj, id, varargin)
         coords = imfrmtIndexToCoordinate(obj.rawDataSize, obj.rawDataFormat, id);
      end
      
      function coords = rawCellIndexToCoordinate(obj, id, varargin)
         coords = imfrmtIndexToCoordinate(obj.rawCellSize, obj.rawCellFormat, id);
      end
      
      

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% casting
  
      function m = maxDataIntensity(obj, varargin)
         m = obj.imaxdataintensity;
      end
      
      function m = minDataIntensity(obj, varargin)
         m = min(obj.data(varargin));
      end
      
      function m = maxRawIntensity(obj, varargin)
         m = obj.imaxrawintensity;
      end
      
      function m = minRawIntensity(obj, varargin)
         m = obj.iminrawintensity;
      end

      function d = castToDataClass(obj, d)
         if iscell(d)
            d = cellfunc(@(x) cast(x,obj.dataclass), d);
         else
            d = cast(d, obj.idataclass);
         end
      end
          
      function d = castToRawClass(obj, d)
         if iscell(d)
            d = cellfunc(@(x) cast(x,obj.rawdataclass), d);
         else
            d = cast(d, obj.irawdataclass);
         end
      end
      
      
      function d = castRawToData(obj, d)
         param.data.min = obj.minRawIntensity;
         param.data.max = obj.maxRawIntensity;
         param.rescale.max = obj.maxDataIntensity;
         param.rescale.min = obj.minDataIntenisty;
         param.class = obj.dataClass;
         
         if iscell(d)
            d = cellfunc(@(x) imrescale(x, param), d);
         else
            d = imrescale(d, param);
         end
      end
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% scaling
      
      function sc = scale(obj)
         dfrmt = 'XYZT';
         pos = imfrmtPosition(dfrmt, obj.dataFormat);
         
         sc = zeros(1,length(pos));
         for i = 1
            sc(i) = obj.iscale(dfrmt(pos(i)));
         end
      end
      
      function sc = scaleX(obj)
         sc = obj.iscale('X');
      end
      
      function sc = scaleY(obj)
         sc = obj.iscale('Y');
      end
 
      function sc = scaleZ(obj)
         sc = obj.iscale('Z');
      end

      function sc = scaleT(obj)
         sc = obj.iscale('T');
      end
      
      function obj = setScale(obj, dim, scale)
         obj.iscale(upper(dim)) = scale;
      end

      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% units
      
      function sc = units(obj)
         dfrmt = 'XYZT';
         pos = imfrmtPosition(dfrmt, obj.dataFormat);
         
         sc = zeros(1,length(pos));
         for i = 1
            sc(i) = obj.iunit(dfrmt(pos(i)));
         end
      end
      
      function sc = unitX(obj)
         sc = obj.iunit('X');
      end
      
      function sc = unitY(obj)
         sc = obj.iunit('Y');
      end
 
      function sc = unitZ(obj)
         sc = obj.iunit('Z');
      end

      function sc = unitT(obj)
         sc = obj.iunit('T');
      end
     
      function obj = setUnit(obj, dim, unit)
         obj.iunit(upper(dim)) = unit;
      end
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% names
      
      function n = name(obj)
         n = obj.iname;
      end
      
      function obj = setName(obj, name)
         obj.iname = name;
      end
          
      function n = channelName(obj, varargin)
         n = obj.ichannelname(varargin{:});
      end
      
      function obj = setChannelName(obj, varargin)
         if nargin == 3
            obj.ichannelname(varargin{1}) = varargin{2};
         else
            obj.ichannelname = varargin{1};
         end
      end
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% meta
      
      function md = metaData(obj, varargin)
         if nargin > 1
            md = obj.imetadata;
            i = find(ismember(md.Parameter, varargin{1}));
            if ~isempty(i)
               md = md.Value(i);
            end
         else
            md = obj.imetadata;
         end
      end 
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% colors
      
      function n = color(obj, varargin)
         n = obj.icolor(varargin{:});
      end
      
      function obj = setColor(obj, varargin)
          if nargin == 3
            obj.icolor(varargin{1}) = varargin{2};
          else
            if ischar(varargin{1})
               obj.icolor = varargin(1);
            else
               obj.icolor = varargin{1};
            end
         end
      end
      
            
      function n = colorAlpha(obj, varargin)
         n = obj.icoloralpha(varargin{:});
      end
      
      function obj = setColorAlpha(obj, varargin)
         if nargin == 3
            obj.icoloralpha(varargin{1}) = varargin{2};
         else
            obj.icoloralpha = varargin{1};
         end
      end

      function obj = reformatColors(obj, ifrmt, ofrmt)
         i = find(lower(ifrmt)== 'c', 1, 'first');
         if ~isempty(i)
            o = find(lower(ofrmt)== 'c', 1, 'first');
            if ~isempty(o)
               if ifrmt(i) ~= ofrmt(o)
                  if ~isempty(obj.icolor)
                     obj.icolor = flip(obj.icolor);
                  end
                  if ~isempty(obj.icoloralpha)
                     obj.ialphacolor = flip(obj.icoloralpha);
                  end
                  if ~isempty(obj.ichannelnames)
                     obj.ichannelnames = flip(obj.ichannelnames);
                  end
               end
            end
         end
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% info

      function istr = infoString(obj, varargin)
         if nargin > 1
            cls = varargin{1};
         else
            cls = '';
         end
         istr = ['ImageSource: ',  cls];
         if ~isempty(obj.name)
            istr = [istr, '\nname:   ', var2char(obj.name)];
         end
         istr = [istr, '\nsize:       ', var2char(obj.dataSize)];
         istr = [istr, '\nformat:     ', var2char(obj.dataFormat)];
         istr = [istr, '\ncellsize    ', var2char(obj.cellSize)]; 
         istr = [istr, '\ncellformat: ', var2char(obj.cellFormat)];     
         istr = [istr, '\nclass:      ', var2char(obj.dataClass)];
         istr = [istr, '\ncolor:      ', var2char(obj.color)];
      end

      function printInfo(obj)
         fprintf([obj.infoString, '\n']);
      end
      
      
   end
end