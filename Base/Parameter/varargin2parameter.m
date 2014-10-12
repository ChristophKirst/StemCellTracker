function param = varargin2parameter(varargin)
%
% param = varargin2parameter(varargin)
%
% description:
%    parses varargin input to a nested parameter struct
%
% input:
%    name, value   pairs, names can contain '.' for nesting
%    param         param struct (may be empty struct)
%    []            
%
% output:
%    param         nested param struct

param = struct();
%fprintf('new call varargin: %s\n', var2char(varargin));
n = 1;

while n <= nargin
   
   in = varargin{n};

   if ~isempty(in) && ~isemptystruct(in)
      
      if ischar(in)
         if nargin < n+1
            error('varargin2parameter: cannot parse parameter, expected value after name: %s !', in);
         end

         %check for sub struct
         pn = strsplit(in, '.');
         param = setsubparam(param, pn, varargin{n+1});
         n = n + 1;
      
      elseif isstruct(in)  
      
         param = appendparameter(param, in);
      
      elseif iscell(in)
      
         addparam = varargin2parameter(in{:});
         param = appendparameter(param, addparam);
      
      else
         if ~isemptyparameter(in)
            error('varargin2parameter: cannot parse parameter, unknonw input: %s !', var2char(in));
         end
      end
   end

   n = n + 1;

end
      

end
      

% helper

function param = setsubparam(param, pn, val)
   pnl = length(pn);
   if pnl == 1
      param.(pn{1}) = val;
      return
   end
   
   if isfield(param, pn{1})
      if ~isstruct(param.(pn{1}))
         error('varargin2parameter: inconsistent sub-stuct structure!')
      end
   else
      param.(pn{1}) = struct();
   end
   
   param.(pn{1}) = setsubparam(param.(pn{1}), pn(2:end), val);
end



function param = appendparameter(param, addparam)
   fns = fieldnames(addparam);
   
   for f = 1:length(fns)
      if ~isfield(param, fns{f})
         param.(fns{f}) = addparam.(fns{f}); %add new field
      else
         if isstruct(param.(fns{f}))
            if ~isstruct(addparam.(fns{f}))
               error('varargin2parameter: inconsistent sub-stuct structure!');
            end
            param.(fns{f}) = appendparameter(param.(fns{f}), addparam.(fns{f}));
         else
            param.(fns{f}) = addparam.(fns{f}); % we dont care if rhs is a struct !
         end          
      end
   end
end


         