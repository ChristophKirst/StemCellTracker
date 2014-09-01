classdef DataSource < matlab.mixin.Copyable
   %
   % DataSource class represents an abstract data interface
   % 
   % decription:
   %     this class can be used to represent data indepenent of its source source
   %     use getData to return the data
   %     use setData to set the data
   %

   properties
      tags = [];   % possible tagging and naming
   end
   
   methods (Abstract)
      dat = getData(obj, varargin);
      suc = setData(obj, varargin);
   end

   methods    
      
      
      
      
      % information      
      function info = getInfo(obj, varargin)
         info.Class = class(obj);
      end
      
      function info = infoString(obj)
         cls = class(obj);
         info = ['DataSource:\nsource: ',  cls(11:end)];          
      end

      function info(obj)
         fprintf([obj.infoString, '\n']);
      end
   end
   
end