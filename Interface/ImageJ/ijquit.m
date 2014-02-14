function ijquit()
%
% mij = ijquit()
%
% description:
%    quits running ImageJ 
%
% See also: ijstart, ijinstance

mij = ijinstance(0);

if ~isempty(mij)
   mij.quit();
end

end

