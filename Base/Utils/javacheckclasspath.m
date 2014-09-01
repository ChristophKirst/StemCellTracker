function check = javacheckclasspath(jarfile, jcpath)
%
% check = javacheckclasspath(name)
%
% description:
%      returns true if jarfile is in the javaclasspath
%
% input:
%      jarfile    the jarfile to check for
%      jcpath     (optional) javaclasspath cell list
%
% output:
%      check       true if jarfile is in classpath
%
% See also: javaddjar, javaclasspath

if nargin < 2
   jcpath = javaclasspath('-all');
end
check = any(cellfun(@(x) ~isempty(regexp(x, ['.*' jarfile '$'], 'once')), jcpath));

end