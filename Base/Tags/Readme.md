Tag Expressions
===============

Replaces tags in a tag expression by numbers and words.
Searches for files with a certain tagformat and returns tag ranges for them
useful to acces (multiple) folders with numbered images 


Example
-------

tf = 'this is <word, s> <count, 2>';
tags.word = 'example';
tags.count = 1;
tagformat2string(tf,  tags)


Version
-------

Version 0.9 
Sep 2014

Copyright:
Christoph Kirst, 
Rockefeller University, 
ckirst@rockefeller.edu

