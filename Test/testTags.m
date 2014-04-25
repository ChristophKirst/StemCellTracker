%%%%%%%%%%%%%%%%%%%%%%%%%
%%% test tagged names %%%
%%%%%%%%%%%%%%%%%%%%%%%%%


%% test some regular expressions
tfrmt = 'test_T<tag1 >-<tag2, 5 >.tif';

tt = regexp(tfrmt, '<(?<name>\w*)\s*(?<k>(,\s*\d*?|\s*?))\s*>', 'names')
{tt.name}
{tt.k}

%%
 tt = regexp(tfrmt, '(?<repl><.*?,*\s*?\d*?\s*?>)', 'names')
 tt.repl

%% tagformat2names
[tagnames, wdt] = tagformat2tagnames(tfrmt)

%% tags2name

name = tags2name(tfrmt, [2, 5])


%% name2tags

name2tags(tfrmt, name)


%% infer tags from stack images

tf = tagformat('./Test/Images/hESCells_tif_movie', {'time','z'})

%%
[tf, tags] = tagformat('./Test/Images/hESCells_tif_movie', {'time','z'})

%%

fns = dir('./Test/Images/hESCells_tif_movie/*');
fns([fns.isdir]) = [];
fns = {fns.name};


regexp(fns, '(?<name>T00\d)\w*?(?<nz>Z0\d)')

%%
name2tags(tf, fns)

%%

tf = tagformat('./Test/Images/hESCells_tif_movie', {'time','z'})

%%
[tf, tags] = tagformat('./Test/Images/hESCells_tif_movie', {'time','z'})

