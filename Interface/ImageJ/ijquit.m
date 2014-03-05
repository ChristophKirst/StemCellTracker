function ijquit()
%
% mij = ijquit()
%
% description:
%    quits running ImageJ 
%
% See also: ijstart, ijinstance

mij = ijinstance();

if ~isempty(mij)
   mij.quit();
end

end

