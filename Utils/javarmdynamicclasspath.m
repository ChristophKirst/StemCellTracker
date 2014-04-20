function javarmdynamicclasspath()
%
% javarmclasspath()
%
% description:
%    tries to remove all dynamic java class paths
%

dyncp = javaclasspath('-dynamic');
javarmpath(dyncp{:});

end