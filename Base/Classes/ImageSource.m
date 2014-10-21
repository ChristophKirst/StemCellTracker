classdef ImageSource < ImageInfo
   %
   % ImageSource class represents abstract Image data 
   % 
   % decription:
   %     this class can be used to represent an Image independent of its actual source
   %     use obj.data to return the actual image data
   %
   % note:
   %     design is as follows: output formats are wrt to the raw formats
   %     a new data / cell data request will read the raw data and reshape it to the output data
   %     if caching is activaed the final cell data will be cached
   %     if raw caching is activated the raw cell data will be cached
   %     for a data request the coordinate ranges are taken over written by specific user cell data specs
   %     the ranges are then inversly reshaed to the raw ranges neccessary to laod the data
    
   properties   
      icache = false;         % (optional) cache the data or not
      icelldata    = [];     % (cached) image data as cell array
      
      irawcache = false;      % (optional) cache the raw cell data
      irawcelldata = [];     % (cached) raw data
      
      idatacorrect = false;      % (optional) switch on/off image correction
      idatacorrectfunction = []; % (optional) function handle to correct images, applied to each data
   end
        
   methods   
      function obj = ImageSource(varargin)  % basic constructor
         %
         % ImageSource()
         % ImageSource(...,fieldname, fieldvalue,...)
         %
         if nargin == 1 
            if isa(varargin{1}, 'ImageSource') %% copy constructor
               obj = copy(varargin{1});
            elseif isnumeric(varargin{1})
               obj = obj.fromData(varargin{1});
            else
               error('%s: not valid arguments for constructor', class(obj));
            end
         else
            obj = obj.fromParameter(varargin);
         end
      end
      
      function obj = fromParameter(obj, varargin)
         obj = classFromParameter(obj, 'i', varargin);
      end
      
      function obj = fromData(obj, data)
         fromData@ImageInfo(obj, data);
         if ~iscell(data)
            data = {data};
         end
         obj.irawcelldata = data; 
      end
   
   
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % cache
      
      function cd = cellDataCache(obj)
         cd = obj.icelldata;
      end
      
      function obj = clearCellDataCache(obj)
         obj.icelldata = [];
      end
      
      function obj = setCellDataCaching(obj, c)
         if c ~= obj.icache
            obj.icache = c;
            obj.clearCellDataCache();
         end
      end

      function cd = rawCellDataCache(obj)
         cd = obj.irawcelldata;
      end
             
      function obj = clearRawCellDataCache(obj)
         obj.irawcelldata = [];
      end
      
      function obj = setRawCellDataCaching(obj, c)
         if c ~= obj.irawcache
            obj.irawcache = c;
            obj.clearRawCellDataCache();
         end
      end
      
            
      function obj = setRawCaching(obj,c)
         % activate raw data caching, usefull for figuring out cell and data formats
         %obj.setCellDataCaching(c);
         obj.setRawCellDataCaching(c); % explicitly specified by user
      end
      
      function obj = setCaching(obj,c)
         % activate caching, usually the final data needs to be cached only
         obj.setCellDataCaching(c);
         %obj.setRawCellDataCaching(c); % explicitly specified by user
      end

      function obj = clearCache(obj)
         obj.clearCellDataCache;
         obj.clearRawCellDataCache;
      end

      function initializeCache(obj)
         % initialize cache by reading all images
         obj.setCaching(true);
         obj.icelldata = obj.cellData; % read all data in to cell
      end
      
      function initializeRawCache(obj)
         % initialize cache by reading all images
         obj.setCaching(true);
         obj.irawcelldata = obj.rawCellData; % read all data in to cell
      end
      

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % data access methods

      % raw data access routine to be overwritten by super class 
      function d = getRawData(obj, rawRange)
         %
         % d = getRawData(obj, varargin)
         %
         % description:
         %     obtain raw image data in the format given by obj.rawformat
         %     range in raw formats and sizes

         cid = obj.rawCellIndexFromRawVarargin(rawRange);
         if length(cid) > 1
            error('%s: getRawData: cell dimensinos are not specified to be singeltons!', class(obj));
         end
         d = obj.irawcelldata{cid}; % trivial here -> place to implement caching
         
         d = imfrmtDataSubset(d, obj.rawDataFormat, rawRange);
      end
 
      % raw cell data access routine to be overwritten by super class 
      function d = getRawCellData(obj, rawRange)
         %
         % d = getRawData(obj, varargin)
         %
         % description:
         %     obtain raw cell data in the format given by obj.rawformat
         %     data specs use raw formats and sizes
         %
         % input:
         %     rawRange    index range restrictions in the raw format and size
         %
         % output:
         %     d           data in rawFormat restricted by raw range
        
                  
         if nargin > 1
            id = obj.rawCellIndexFromRawVarargin(rawRange);
            d  = obj.irawcelldata(id);     
         else
            d = obj.irawcelldata; % trivial here -> place to implement caching
         end
         
         % restrict the data for range
         for i = 1:length(d)
            d{i} = imfrmtDataSubset(d{i}, obj.rawDataFormat, rawRange);
         end
      end
 

      function d = getDataFromRaw(obj, varargin)  % obtain the image data
         %
         % d = getData(obj, varargin)
         %
         % description:
         %     obtain data from raw data
         %     determines inverse coordinate range
         %     then uses getRawCellData to obtain raw data 
         %     then reshape to output data
         %     checks result for single data cell
         %
         % note: indices passed as arguments refer to cell indices, and full data of that cell is returned
        
         % check if all cell dims are singeltons
         cid = obj.cellIndex(varargin{:});
         if length(cid) > 1
            error('%s: getData: cell dimensions not specified to singeltons!', class(obj));
         end
         
         
         % ranges necessary to load the raw data
         [rawRange, rawReshapeSize] = obj.rawRange(varargin{:});
         
         % load raw data as cell data
         cd = obj.getRawCellData(rawRange);
         
           
         % transform to output cell data 
         fullDataFrmt = obj.fullDataFormat;
         fullCellFrmt = obj.fullCellFormat;
         
         cd = imfrmtReshapeCellData(cd, obj.rawDataFormat, obj.rawCellFormat, fullDataFrmt, fullCellFrmt, ...
                                        obj.reshapeFrom, obj.reshapeTo, rawReshapeSize);
                                     
         cd = imfrmtReformatCellData(cd, fullDataFrmt, fullCellFrmt, obj.dataFormat,  obj.cellFormat);

         % cd should be singeston
         if length(cd) ~= 1
            error('%s: getDataFromRaw: incsonsistent data sizes!', class(obj));
         end
         
         d = cd{1};
      end
      
      
      function cd = getCellDataFromRaw(obj, varargin)
         %
         % cd = getCellDataFromRaw(obj, varargin)
         %
         % description:
         %     determines inverse coordinate range
         %     then uses getRawCellData to obtain raw cell data 
         %     then reshape to output cell data

         % ranges necessary to load the raw data
         [rawRange, rawReshapeSize] = obj.rawRange(varargin{:});
         
         % load raw data as cell data
         cd = obj.getRawCellData(rawRange);
         
         % transform to output cell data 
         fullDataFrmt = obj.fullDataFormat;
         fullCellFrmt = obj.fullCellFormat;
         
         cd = imfrmtReshapeCellData(cd, obj.rawDataFormat, obj.rawCellFormat, fullDataFrmt, fullCellFrmt, ...
                                        obj.reshapeFrom, obj.reshapeTo, rawReshapeSize);
                                     
         cd = imfrmtReformatCellData(cd, fullDataFrmt, fullCellFrmt, obj.dataFormat,  obj.cellFormat);
         
      end
  
      function obj = setData(obj, d, varargin)
         id = obj.dataIndex(varargin{:});
         if length(id) ~= 1
            error('%s: setData: cell dimensions not specified to singeltons!', class(obj));
         end 
         obj.icelldata{id} = d;
      end
      
      
      function obj = setRawData(obj, d, varargin)
         id = obj.rawDataIndex(varargin{:});
         if length(id) ~= 1
            error('%s: setRawData: cell dimensions not specified to singeltons!', class(obj));
         end 
         obj.irawcelldata{id} = d;
      end

      
      function obj = setCellData(obj, cd)
         obj.icelldata = cd;
      end
      
      function obj = setRawCellData(obj, cd)  % set the image data
         obj.irawcelldata = cd;
      end
      

      
      
 
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % data methods      

      function d = data(obj, varargin)
         if obj.icache % basic caching
            
            id = obj.cellIndex(varargin{:});
            if length(id) ~= 1
               error('%s: data: cell size not a singelton!');
            end
            
            if isempty(obj.icelldata)
               obj.icelldata = cell(obj.cellSize);
            end
               
            if ~isempty(obj.icelldata{id})
               d = obj.icelldata{id};
            else
               d = obj.getDataFromRaw(varargin{:});
               
               % image correction 
               d = obj.correctData(d);
               
               obj.icelldata{id} = d;
            end

         else
            d = obj.getDataFromRaw(varargin{:});
            
            % image correction 
            d = obj.correctData(d);
         end
      end
      
      
      function c = cell(obj, varargin)
         if obj.icache % basic caching
            
            id = obj.cellIndex(varargin{:});

            if isempty(obj.icelldata)
               obj.icelldata = cell(obj.cellSize);
            end

            idload = cellfun(@isempty, obj.icelldata(id));
            idload = id(idload);
               
            for i = 1:length(idload)
               obj.icelldata{idload(i)} = obj.data(idload(i));
            end
            
            c = obj.icelldata(id);
            
            
         else
            c = obj.getCellDataFromRaw(varargin{:});
            
            % (optional) image correction
            c = cellfunc(@(x)obj.correctData(x), c);
         end
      end
      
          
      function d = dataSubset(obj, varargin)
         %
         % d = dataSubset(obj, datasepc)
         %
         % description:
         %    subset of the data given the data specifications datasepc
         %
         
         d = imfrmtDataSubset(obj.data, obj.dataFormat, varargin);
      end
      
      function d = dataExtract(obj, roi)
         %
         % d = subdata(obj, datasepc)
         %
         % description:
         %    extract a subset of the data given the spatial roi 
         %
         
         d = roi.extractData(obj.data); 
      end

      function d = dataResample(obj, scalefac, varargin)
         d = obj.data(varargin{:});
         d = imresize(d, scalefac);
      end
      
      function d = cellResample(obj, scalefac, varargin)
         d = obj.cell(varargin{:});
         d = cellfun(@(x) imresize(x, scalefac), d);
      end

  
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% reshaping / reformatting data     
      
%       function obj = setDataFormat(obj, dfrmt)         
%          if ~isempty(obj.icelldata)
%             obj.icelldata = obj.reformatData(obj.icelldata, obj.dataFormat, dfrmt);
%          end
%          setDataFormat@ImageInfo(obj, dfrmt);
%       end
      
 
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% casting and intensity scaling
     
      function obj = initializeRawDataIntensity(obj, varargin)
         cd = obj.rawCellData(varargin{:});
         m = cellfun(@(x) min(x(:)), cd);
         obj.iminrawintensity = min(m(:));
         m = cellfun(@(x) max(x(:)), cd);
         obj.imaxrawintensity = max(m(:));  
      end
      
      function obj = initializeDataIntensity(obj, varargin)
         cd = obj.cell(varargin{:});
         m = cellfun(@(x) min(x(:)), cd);
         obj.iminintensity = min(m(:));
         m = cellfun(@(x) max(x(:)), cd);
         obj.imaxintensity = max(m(:));   
      end 
      
          
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% image correction
     
      function obj = setDataCorrect(obj, c)
         obj.idatacorrect = c;
      end
      
      function obj = resetDataCorrect(obj)
         obj.idatacorrect = false;
         obj.idatacorrectfunction = [];
      end
 
      function c = dataCorrect(obj)
         c = obj.idatacorrect;
      end

      function c = dataCorrectFunction(obj)
         c = obj.idatacorrectfunction;
      end
      
      function obj = setDataCorrectFunction(obj, c)
         obj.idatacorrectfunction = c;
      end
      
      function obj = setBackgroundAndFlatFieldCorrection(obj, bkg, flt)
         if nargin < 3
            flt = ones(sizebkg);
         end
         obj.idatacorrectfunction = @(x) correctFromBackgroudAndFlatField(x, bkg, flt);
         obj.idatacorrect = true;
      end
      
      function d = correctData(obj, d)
         if obj.dataCorrect
            c = obj.idatacorrectfunction;
            if ~isempty(c)
               d = c(d);
            end
         end
      end
      
      function [bkg, flt] = backgroundAndFlatFieldFromMin(obj, varargin)
         [bkg, flt] = backgroundFromMin(obj, varargin{:});
      end

      function [bkg, flt] = backgroudnAndFlatFieldFromDataCorrectFunction(obj)
         [bkg, flt] = backgroundAndFlatFieldFromCorrectBackgroundAndFlatFieldHandle(obj.idatacorrectfunction);
      end
      
      function obj = initializeBackgroundCorrectionFromMin(obj, varargin)
         [bkg, flt] = backgroundFromMin(obj, varargin{:});
         obj.idatacorrectfunction = @(x) correctFromBackgroudAndFlatField(x, bkg, flt);
      end
      
      function obj = initializeBackgroundCorrectionFromBackgroudnAndFlatField(obj, bkg, flt)
         obj.idatacorrectfunction = @(x) correctFromBackgroudAndFlatField(x, bkg, flt);
      end
      

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% alignment utils
 
      function d = cellPreview(obj, varargin)
         
         d2 = alignCellPreview(obj, varargin);

         if nargout < 1
            implot(d2)
         else
            d = d2;
         end
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % information / visulaization 
      
      function plot(obj)
         plotImageSource(obj);
      end
      
      function plottiling(obj, varargin)
         implottiling(obj.cell(varargin{:}), varargin)
      end
      
   end
   
end