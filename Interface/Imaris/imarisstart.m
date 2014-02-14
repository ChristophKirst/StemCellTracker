function imaris = imarisstart(varargin)
%
% imaris = istart(ipath)
% imaris = istart(id)
% imaris = istart(ipath, id);
%
% description:
%    starts Imaris application if not already running
%
% input:
%   ipath   (optional) path to Imaris (autodetect)
%   id      (optional) id of Imaris application (0)
%
% output:
%   imaris  reference to Imaris application
%
% See also: imarisinitialize, imarisid, isimarisid

id = {};
ipath = {};

if nargin > 0
   if isnumeric(varargin{1})
      id = varargin(1);
   elseif ischar(varargin{1})
      ipath = varargin(1);
      if nargin > 1
         id = varargin(2);
      end
   end
end

% check for running instance
imaris = imarisinstance(id{:});
if ~isempty(imaris)
   return
end
id = 42;

[ipath, server, exe] = imarispath(ipath{:});
if isempty(ipath)
   error('imarisstart: cannot find Imaris.exe')
end

startsever(server);

imaris = startimaris(exe, id);

end



%%% Helper

function imaris = startimaris(exe, id)

[status, result] = system(['"', exe, '" id', num2str(id), ' &']);
if status == 1
   error('imarisstart: error starting Imaris: %s', result)
end

% give it some time
wait(1);

try
   for trial = 1 : 200
      try
         ilib = ImarisLib();
         imaris = ilib.GetApplication(id);
      catch ex %#ok<NASGU>
         % We can't catch the Imaris exception...
      end
      
      if ~isempty(imaris)
         break
      end
      pause(0.1)
   end
   
   if isempty(imaris)
      error('imarisstart: error linking to Imaris')
   end
   
catch ex
   error('imarisstart: error: %s', ex.message)
end

end



function startsever(server)

   % check for running ImarisServerIce.exe
   if isserverrunning()
      return
   end
   
   % start instance and wait until it is running
   try
      % Launch ImarisServer
      [status, result] = system(['"', server, '" &']);
    
      if status == 1
         error('imarisstart: error when launching server:\n%s', result);
      end
    
      % wait for server
      tic;
      t = 0;
      timeout = 10;
      while t < timeout
        
         if isserverrunning() == 1
            return;
         end
        
         t = toc;  
      end
      
      error('imarisstart: timeout for starting the Imaris server\n%s', server);
    
   catch ex
      
      error('imarisstart: error: %s', ex.message);
   end

end


function b = isserverrunning()
   b = 0;
   if ispc()
      [~, result] = system(...
         'tasklist /NH /FI "IMAGENAME eq ImarisServerIce.exe"');
      if strfind(result, 'ImarisServerIce.exe')
         b = 1;
         return;
      end

   elseif ismac()
      [~, result] = system('ps aux | grep ImarisServerIce');
      if strfind(result, this.mImarisServerExePath)
         b = 1;
         return;
      end

   else
      error('imarisstart: unsupported platform.');
   end
end

