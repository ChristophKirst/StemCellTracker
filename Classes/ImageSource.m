classdef ImageSource < matlab.mixin.Copyable
   %
   % ImageSource class represents abstract Image data 
   % 
   % decription:
   %     this class can be used to represent an Image independent of its actual source
   %     use getData to return the actual image data
   %
   properties
      size   = [];                  % size of the image
      format = '';                  % pql format of the image
      
      colors = [];                  % colors of the channels
      tagnames = containers.Map();  % translation between names and tags, e.g. for channels 'dapi' -> ('channel' -> 1)
   end
     
   methods (Abstract)
      img = getRawData(obj, varargin);
      siz = getSize(obj, varargin);
   end
   
   
   methods    
      function initialize(obj)
         if isempty(obj.size)
            obj.size = obj.getSize();
         end
         if isempty(obj.format)
            obj.format = obj.getFormat();
         end
      end

      function frmt = getFormat(obj, varargin) 
         frmt = imsize2format(obj.size);
      end
      
      function info = getInfo(obj, varargin)
         info.Class = class(obj);
         info.Format = obj.format;
         info.Size = obj.size;
      end
      
      function img = getData(obj, varargin)
         img = obj.getRawData(varargin{:});
         img = impqlpermute(img, obj.format, 'pqlct');
      end
      
      
      
      
      % information
      function info = infoString(obj)
         cls = class(obj);
         info = ['ImageSource:\nsource: ',  cls(12:end), '\nsize:   ', var2char(obj.size), '\nformat: ', var2char(obj.format)]; 
         
         if ~isempty(obj.tagnames)
            info = [info '\ntagnames:'];
            keys = obj.tagnames.keys;
            vals = obj.tagnames.values;
            for i = 1:obj.tagnames.length
               info = [info, var2char(keys{i}), ' -> ' var2char(vals{i})]; %#ok<AGROW>
            end
         end
         
         
         
      end

      function info(obj)
         fprintf([obj.infoString, '\n']);
      end
   end
   
end