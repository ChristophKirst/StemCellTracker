function illoadclassifier(clfile)
%
% illoadclassifier(clfile)
%
% description:
%    set up the Ilastik clasifier by loading it form the h5 file clfile
%
% input:
%    clfile   classifier file as generated with ilastic software
%
% See also: ilinitialize, ilclassify

py('eval', ['ilc.load_classifier(''' clfile ''')'])



