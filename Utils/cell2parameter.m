function param = cell2parameter(varargin)
%
% param = cell2parameter(varargin)
% param = cell2parameter(param, varargin)
%
% description
%    returns a struct with parameters given by varargin
%    paramter names 'xxx.yyy' are parased to sub structs
%
% input:
%   varargin    parameter pairs parname, parval
%   param       existing parameter struct 
% 
% ouput:
%  param        parameter struct


if nargin == 0
   param = [];
   return
end

if isstruct(varargin{1})
   param = varargin{1};
   varargin = varargin(2:end);
   narg = nargin - 1;
else
   param = struct();
   narg = nargin;
end

if mod(narg, 2) ~= 0
   error('setParameter: expect even number of arguments!');
end

for i = 1:2:narg
   pn = varargin{i};
   
   %check for sub struct
   if ischar(pn)
      pn = strsplit(pn, '.');
   elseif ~iscellstr(pn)
      error('setParameter: expect parameter name or cellstr for input at pos %g!', i);
   end
   param = setsubparam(param, pn, varargin{i+1});
end

end
      

function par = setsubparam(param, pn, val)
   pnl = length(pn);
   par = param;
   if pnl == 1
      par.(pn{1}) = val;
      return
   end

   if isfield(par, pn{1})
      if ~isstruct(par.(pn{1}))
         error('setParameter: inconsistent sub stuct structure!')
      end
   else
      par.(pn{1}) = struct();
   end
   
   par.(pn{1}) = setsubparam(par.(pn{1}), pn(2:end), val); 
end


         
      


