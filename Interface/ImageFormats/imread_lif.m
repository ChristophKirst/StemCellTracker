function data = imread_lif(name, varargin)
%
% data = imread_lif(name, id)
%
% description: 
%     reads imagedata from lif file using bio-formats tools
%
% input:
%     name   filename
%     id     (optional) series id, if not given the gui is opened
%
% output"
%     data   image data arrray
%
% See also: bfopen

data = imread_bf(name, varargin{:});

end

