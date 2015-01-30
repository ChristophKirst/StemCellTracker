function ijrun(cmd, varargin)
%
% ijrun(cmd, varargin)
%
% description:
%    Runs a ImageJ macro command cmd with posible option varargin
%
% input:
%    cmd       command to run as string
%    varargin  options to command as string

if nargin > 1
   MImageJ.run(cmd, varargin{1});
else
   MImageJ.run(cmd);
end

end
   