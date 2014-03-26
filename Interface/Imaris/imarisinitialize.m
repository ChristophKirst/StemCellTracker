function imarisinitialize(ipath)
%
% imarisinitialize(ipath)
%
% description:
%    initializes Imaris enviroment by making ImarisLib.jar from Imaris XT
%    visible to matlab
%
% input:
%    ipath  (optional) path to Imaris ([] = automatic detection)
%
% See also: imarisstart, imarisinstance

% Imaris only runs on Windows and Mac OS X
if ~(ismac() || ispc())
    error('imarisinitialize: Imaris runs only on Windows or Mac OS X.');
end

if nargin < 1 || isempty(ipath)
   ipath = {};
else
   ipath = {ipath};
end


% Check if Imaris java libs are already installed 
jcp = javaclasspath('-all');
jpc = strfind(jcp, 'ImarisLib.jar');
jpc = ~cellfun(@isempty, jpc);

if any(jpc)  % ImarisLib.jar already installed
   fprintf('imarisinitialize: ImarisLib.jar already installed.\n'); 
   %success = true;
   return;
end

% Search for ImarisLib.jar
ipath = imarispath(ipath{:});

if isempty(ipath)
   error('imarisinitialize: could not find ImarisLib.jar.'); 
   %return;
else
   %add ImarisLib.jar to java libraries
   libjar = fullfile(ipath, 'XT', 'matlab', 'ImarisLib.jar');
   success = javaaddjar(libjar);
   if success; 
      fprintf('imarisinitialize: installed %s\n', libjar); 
   else
      error('imarisinitialize: could not find ImarisLib.jar.'); 
      %return;
   end
end


% On Mac OS X, make sure that the Imaris Frameworks folder is at the
% beginning of the dynamic library path to prevent any conflicts with
% the Qt libraries loaded by MATLAB
if ismac()
    dylPath = getenv('DYLD_LIBRARY_PATH');
    imarisFrameworksPath = [ipath, filesep, 'Contents', filesep, 'Frameworks'];
    indx = strfind(dylPath, imarisFrameworksPath);
    if isempty(indx) || indx ~= 1
        % Prepend
        dylPath = [imarisFrameworksPath, ':', dylPath];
        setenv('DYLD_LIBRARY_PATH', dylPath)
    end
end


end