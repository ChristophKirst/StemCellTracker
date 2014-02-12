function imaris = imarisstart(ipath)
%
% imaris = istart(ipath)
%
% description:
%    starts Imaris searching for XT in ipath
%
% input:
%    ipath   path to ImarisXT or id of instance
%
% notes:
%    based on IceImarisConnector


% Imaris only runs on Windows and Mac OS X
if ~(ismac() || ispc())
    error('imarisstart: Imaris can only work on Windows and Mac OS X.');
end

if nargin < 1
   ipath = 0;
end

% On Mac OS X, make sure that the Imaris Frameworks folder is at the
% beginning of the dynamic library path to prevent any conflicts with
% the Qt libraries loaded by MATLAB
if ismac()
    dylPath = getenv('DYLD_LIBRARY_PATH');
    ImarisFrameworksPath = [this.mImarisPath, filesep, ...
        'Contents', filesep, 'Frameworks'];
    indx = strfind(dylPath, ImarisFrameworksPath);
    if isempty(indx) || indx ~= 1
        % Prepend
        dylPath = [ImarisFrameworksPath, ':', dylPath];
        setenv('DYLD_LIBRARY_PATH', dylPath)
    end
end






if isscalar(ipath)
  
                   
     % Check if the application is registered
     server = mImarisLib.GetServer();
     if isempty(server)
         error('Could not connect to Imaris Server.');
     end
     nApps = server.GetNumberOfObjects();
     if nApps == 0
         error('There are no registered Imaris applications.');
     end

     % Does the passed ID match the ID of any of the
     % registered (running) Imaris application?
     found = 0;
     for i = 0 : nApps - 1
         if server.GetObjectID(i) == imarisApplication
             found = 1;
             break;
         end
     end
     if found == 0
         error('Invalid Imaris application ID.');
     end

     % Now we can get the application object and store it
     imarisApplicationObj = ...
         this.mImarisLib.GetApplication(imarisApplication);
     if isempty(imarisApplicationObj)
         error('Cannot access Imaris application.');
     else
         this.mImarisApplication = imarisApplicationObj;
     end

   else

     error(['The passed object is not a Imaris ', ...
         'Application ID.']);

   end
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   ipath = findImarisPath();
end















% Now we open a new one
try
    
    % Start ImarisServer
    [success, errorMessage] = this.startImarisServer();
    if success == 0
        disp(errorMessage);
        success = 0;
        return
    end
    
    % Launch Imaris
    [status, result] = system(...
        ['"', this.mImarisExePath, '" id', ...
        num2str(this.mImarisObjectID), ' &']);
    if status == 1
        disp(result);
        success = 0;
        return
    end
    
    % Try getting the application over a certain time period in case it
    % takes to long for Imaris to be registered.
    for trial = 1 : 200
        try
            % A too quick call to mImarisLib.GetApplication() could potentially
            % throw an exception and leave the mImarisLib object in an unusable
            % state. As a workaround, we reinstantiate ImarisLib() at 
            % every iteration. This will make sure that sooner or later 
            % we will get the application. The exception is not
            % automatically casted to a MATLAB exception, se we cannot 
            % really catch it...
            this.mImarisLib = ImarisLib();
            vImaris = this.mImarisLib.GetApplication(this.mImarisObjectID);
        catch ex %#ok<NASGU>
            % We can't catch the exception...
        end
        if ~isempty(vImaris) 
            break 
        end
        pause(0.1)
    end
    
    % At this point we should have the application
    if isempty(vImaris) 
        disp('Could not link to the Imaris application.');
        success = 0;
        return
    end
    
    % We can store the application
    this.mImarisApplication = vImaris;

    % We set success to 1
    success = 1;
    
catch ex
    
    disp(ex.message);
    success = 0;
    
end

end












if ischar(ipath)
   
   
   
   
   
elseif isnumeric(ipath)
   
else
   error('imarisstart: ipath should be char or instance id')
end
   
   



         % First, we prepare everything we need
            % Store the Imaris and ImarisLib path
            [success, errorMessage] = this.findImaris;
            if ~success
                error(errorMessage);
            end
            
            % Add the ImarisLib.jar package to the java class path
            % (if not there yet)
            if all(cellfun(...
                    @isempty, strfind(javaclasspath, 'ImarisLib.jar')))
                javaaddpath(this.mImarisLibPath);
            end
            
            % Create and store an ImarisLib instance
            this.mImarisLib = ImarisLib();

            % Assign a random id
            this.mImarisObjectID = randi(100000);

            % Now we check the (optional) input parameter.
            % If the constructor is called without parameters, we just
            % create an IceImarisConnector object that does nothing.
            % If we get one input parameter, we have to distinguish between
            % three cases: 
            % - we get an Imaris Application ID as provided by Imaris and
            %      thus we query the Imaris Server for the application
            % - we get an IceImarisConnector reference and thus we just
            %      return it
            % - we get an Imaris Application ICE object (rare) and thus
            %      we simply assign it to the mImarisApplication property.
            
            if nargin == 0
                
                % We already did everything
                return

            end
                
            if nargin > 0 && ~isempty(imarisApplication)
                
                if isa(imarisApplication, 'IceImarisConnector')
                    
                    % If the input parameter is an IceImarisConnector
                    % object we return the reference. This way, an
                    % XTension class can take a reference to an 
                    % IceImarisConnector object as input parameter
                    this = imarisApplication;
                    
                elseif isa(imarisApplication, ...
                        'Imaris.IApplicationPrxHelper')
                    
                    % This is an Imaris application object - we store it
                    this.mImarisApplication = imarisApplication;
                    
    
                
            end
            
            if nargin > 1 && ~isempty(indexingStart)
                
                this.mIndexingStart = indexingStart;
                
            end
            
            if nargin > 2
                
                error('Wrong number of input arguments!');
                
            end
            
        end























%% Subfunctions

function [ipath, status, errorMessage] = findImarisPath(this)

status = 0;
errorMessage = '';

ipath.mImarisPath = '';
ipath.mImarisExePath = '';
ipath.this.mImarisServerExePath = '';
this.mImarisLibPath = '';


% Try to get the path from the environment variable IMARISPATH
imarisPath = getenv('IMARISPATH');



% Is the variable defined?
if isempty(imarisPath)
    % jonas - 4/2012 be a little more robust
    % installation defaults to C: on Windows, /Applications on Mac
    if ispc
        tmp = 'C:\Program Files\Bitplane\';
    elseif ismac
        tmp = '/Applications';
    else
        errorMessage = ...
            'IceImarisConnector can only work on Windows and Mac OS X';
        return
    end
    
    if exist(tmp,'dir')
        % Pick the directory name with highest version number
        % Aaron - 09/12. Use highest version number instead of
        %                latest modification date to choose.
        d = dir(fullfile(tmp,'Imaris*'));
        newestVersionDir = findNewestVersion(d);
        if isempty(newestVersionDir)
            errorMessage = sprintf(...
                ['No Imaris installation found in %s.',...
                ' Please define an environment variable ', ...
                'IMARISPATH'],...
                tmp);
            return;
        else
            imarisPath = fullfile(tmp, newestVersionDir);
        end
    else
        errorMessage = sprintf(...
            ['No Imaris installation found in %s.',...
            ' Please define an environment variable IMARISPATH'],...
            tmp);
        return;
    end
    
else

    % Does it point to an existing directory?
    if ~exist(imarisPath, 'dir')
        errorMessage = ['The content of the IMARISPATH environment ', ...
            'variable does not point to a valid directory.'];
        return;
   end
    
end

% Now store imarisPath and proceed with setting all required executables
% and libraries
this.mImarisPath = imarisPath;

% Set the path to the Imaris and ImarisServer executables, and to the 
% ImarisLib library
if ispc()
    exePath = fullfile(imarisPath, 'Imaris.exe');
    serverExePath = fullfile(imarisPath, 'ImarisServerIce.exe');
    libPath = fullfile(imarisPath, 'XT', 'matlab', 'ImarisLib.jar');
elseif ismac()
    exePath = fullfile(imarisPath, 'Contents', 'MacOS', 'Imaris');
    serverExePath = fullfile(imarisPath, 'Contents', 'MacOS', 'ImarisServerIce');
    libPath = fullfile(imarisPath, 'Contents', 'SharedSupport', ...
        'XT', 'matlab', 'ImarisLib.jar');
else
    errorMessage = ['IceImarisConnector can only be used on Windows ', ...
        'and Mac OS X.'];
    return
end

% Check whether the executable Imaris file exists
if ~exist(exePath, 'file')
    errorMessage = 'Could not find the Imaris executable.';
    return;
end

% Check whether the executable ImarisServer file exists
if ~exist(serverExePath, 'file')
    errorMessage = 'Could not find the ImarisServer executable.';
    return;
end

% Check whether the ImarisLib jar package exists
if ~exist(libPath, 'file')
    errorMessage = 'Could not find the ImarisLib jar file.';
    return;
end

% Now we can store the information and return success
this.mImarisExePath = exePath;
this.mImarisServerExePath = serverExePath;
this.mImarisLibPath = libPath;

status = 1;

% In case of multiple Imaris installations, return the most recent
function newestVersionDir = findNewestVersion(allDirs)
        
        % If found, this will be the (relative) ImarisPath
        newestVersionDir = [];
        
        % Newest version. Initially set to one since valid versions will
        % be larger, invalid versions might be zero.
        newestVersion = 1;
        
        % Make sure to ignore the Scene Viewer, the File Converter and 
        % the 32bit version on 64 bit machines
        allDirs(~cellfun(@isempty, ....
            strfind({allDirs.name}, 'ImarisSceneViewer'))) = [];
        allDirs(~cellfun(@isempty, ....
            strfind({allDirs.name}, 'FileConverter'))) = [];
        allDirs(~cellfun(@isempty, ....
            strfind({allDirs.name}, '32bit'))) = [];
        
        for i = 1 : numel(allDirs)
            
            % Extract version from directory name
            tokens = regexp(allDirs(i).name, ...
                '(\d)+\.(\d)+\.+(\d)?', 'tokens');
            
            if isempty(tokens) || numel(tokens{1}) ~= 3
                continue;
            end
            
            % Get the major, minor and patch versions
            major = str2double(tokens{1}{1});
            if isnan(major)
                % Must be defined
                continue;
            end
            
            minor = str2double(tokens{1}{2});
            if isnan(minor)
                % Must be defined
                continue;
            end
            
            patch = str2double(tokens{1}{3});
            if isnan(patch)
                % In case the patch version is not set we assume 0 is meant
                patch = 0;
            end
            
            % Compute version as integer
            version = 1e6 * major + 1e4 * minor + 1e2 * patch;
            
            % Is it the newest yet?
            if version > newestVersion
                newestVersionDir = allDirs(i).name;
                newestVersion = version;
            end
        end

    end
end

