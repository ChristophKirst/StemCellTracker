function p = initializeParallelProcessing(numWorkers)
%
% intializes the parallel 
%

if nargin < 1
    numWorkers = {};
else
    numWorkers = {numWorkers};
end

if isempty(gcp('nocreate'))
   p = parpool(numWorkers{:});
    
else
   p = gcp;
   if ~isempty(numWorkers) && numWorkers{1} ~= p.NumWorkers
      delete(p)
      p = parpool(numWorkers{:});
   end
end

% set java class path in the workers
clp = javaclasspath;
spmd
    initialize
    javaclasspath(clp)
end

