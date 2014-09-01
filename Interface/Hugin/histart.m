function histart(varargin)
%
% histart()
%
% description:
%    starts hugin

if nargin > 0 && ischar(varargin{1})
   cmd = ['hugin ' varargin{1}, ' &'];
else
   cmd = 'hugin &';
end

system(cmd)

end