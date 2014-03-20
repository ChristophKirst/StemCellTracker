function v = bfversion()
%
% v = bfversion()
%
% description:
%    returns version of bf tools (loci_tools.jar)
%
% output:
%    v   version string
%
% See also: bfinitialize

v = loci.formats.FormatTools.VERSION;
v = char(v);

end