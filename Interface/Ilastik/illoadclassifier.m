function ilc = illoadclassifier(clfile)
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

ilc = py.IlastikInterface.IlastikClassifier();
ilc.load_classifier(clfile);

end


