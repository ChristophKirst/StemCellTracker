function [bkg, flt] = backgroundAndFlatFieldFromCorrectBackgroundAndFlatFieldHandle(h)
%
% [bkg, flt] = backgroundAndFlatFieldFromCorrectBackgroundAndFlatFieldHandle(h)
%
% description:
%     restores backgroudn anf flatfield form a function handle of the form
%     @(x)(correctFromBackgroudAndFlatField(x,bkg,flt)


f = functions(h);

% idellay: parse the bkg and flatfield exprs here and eval them using the workspace 
%fn = f.function;

bkg = f.workspace{1}.bkg;
flt = f.workspace{1}.flt;

end