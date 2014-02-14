Imaris Interface
================

Interface to Imaris to visualize anaylsis results


Usage
-----

To start imaris run

   imarisinitialize
   imarisstart

open file

   imarisopen(filename)

see whats there

   objects = imarisobjects()

get/set data

   volume = imarisget(objectname)                     % volumetirc data
   [vertices, faces, normals] = imarisget(objectname) % surface data

   imarisset(volumedata)                              % will find the first availabel volume
   imarisset(objectname, vertices, faces, normals)

save

   imarissave(filename)

quit

   imarisquit();


for full functionality see the individual m files


Author
------

Christoph Kirst
The Rockefeller University, New York, 2014


Acknowledgements
----------------

inspired by IceImarisConnector


