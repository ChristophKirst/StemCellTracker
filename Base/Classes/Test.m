classdef Test < handle
   
   properties
      x = struct('v', 1, 'w', 2, 'z', 60);
      y = [];
      
      
   end
   
   methods 
      function obj = Test(varargin) %xx,yy)
         if nargin > 0
            properties(obj)
         end
      end
      
      
      function d = v(obj, varargin)
         nargin
         
         %d = [obj.x];
         %d = [d.v];
      end
         
      
   end
end
      