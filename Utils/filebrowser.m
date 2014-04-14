function filebrowser(filename)
%
% filebrowser(filename)
%
% description:
%     opens filename in systems default file browser
%
% input:
%     filename   directory or filename
%

if nargin < 1
   filename = '.';
end

fdir = fileparts(filename);
if isempty(fdir)
   fdir  = '.';
end

% Windows PC    
if ispc
    C = evalc(['!explorer ' fdir]);

% Unix or derivative
elseif isunix

    % Mac
    if ismac
        C = evalc(['!open ' fdir]);

    % Linux
    else
        fMs = {...
            'xdg-open'   % most generic one
            'gvfs-open'  % successor of gnome-open
            'gnome-open' % older gnome-based systems               
            'kde-open'   % older KDE systems
           };
        C = '.';
        ii = 1;
        while ~isempty(C)                
            C = evalc(['!' fMs{ii} ' ' fdir ' &']);
            ii = ii +1;
        end

    end
else
    error('filebrowser: Unrecognized operating system.');
end

if ~isempty(C)
    error(['filebrowser: Error while opening directory in default file manager.\n',...
        'The reported error was:\n%s'], C); 
end