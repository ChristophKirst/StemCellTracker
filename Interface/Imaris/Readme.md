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

   volume = imarisgetvolume()                                % volumetirc data
   [vertices, faces, normals] = imarisgetsurface(objectname) % surface data

   imarissetvolume(volumedata)                               
   imarissetsurface(objectname, vertices, faces, normals)

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


