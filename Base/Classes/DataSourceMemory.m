classdef DataSourceMemory < DataSource
   %
   % DataSourceMemory represents data in memory
   % 
 
   properties 
      data = [];  % the data array
   end

   methods
      
      function dat = getData(obj, varargin)
         dat = obj.data;
      end
      
      function suc = setData(obj, varargin)
         obj.data = varargin{1};
         suc = true;
      end
      
   end
  
end