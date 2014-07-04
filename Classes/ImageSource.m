classdef ImageSource < matlab.mixin.Copyable
   %
   % ImageSource class represents abstract Image data 
   % 
   % decription:
   %     this class can be used to represent an Image independent of its actual source
   %     use readImage to return the actual image data
   %
   properties
      size   = [];
      format = '';
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
      
      function info = infoString(obj)
         cls = class(obj);
         info = ['ImageSource:\nsource: ',  cls(12:end), '\nsize:   ', var2char(obj.size), '\nformat: ', var2char(obj.format)]; 
      end

      function info(obj)
         fprintf([obj.infoString, '\n']);
      end
   end
   
end