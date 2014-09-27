Parameter Handling
==================

introduces easy to read and code systematic parameter access / passing / parsing 


Example
-------

function r = y(x, varargin)
   param = parseParameter(varargin);
   r = getParameter(param, 'return', 'hello world');
end



Version
-------

Version 0.9 
Sep 2014

Copyright:
Christoph Kirst, Rockefeller University, 
ckirst@rockefeller.edu
