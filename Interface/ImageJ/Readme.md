ImageJ Interface
================

Provides a link to ImageJ

Features:
---------
   -ijplot3d    plot 3d graphics
   -ijget       get current image openin imagej
   -ijread      read images via imagej libraries
   -ijimage2mat converts ImagePlus to matlab array

   
Requirements:
-------------
   - A working ImageJ or Fiji installation
   - java's j3d


Installation:
-------------
   - to compile java class MImageJ.java run ijcompile.m
   - to install java's j3d use the InstallJava3D script in /Utils/External/Fiji


Note:
-----

Matlab generates errors of the form 
Undefined variable ... or ... if initializaion of a class fails due tio missing 
java libraries, such as j3d !