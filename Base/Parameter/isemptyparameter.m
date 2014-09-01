function ie = isemptyparameter(param)
%
% ie = isemptyparameter(param)
%
% description:
%    checks if param isempty or empty struct
%
% input:
%    param   parameter struct or []
%
% output:
%   ie       true if param is [] or empty struct
%
% note: matlabs isempty returns 1 for struct()
%
% See also: isempty

if isstruct(param) &&  isempty(fieldnames(param))
   ie = true;
else
   ie = isempty(param);
end

end