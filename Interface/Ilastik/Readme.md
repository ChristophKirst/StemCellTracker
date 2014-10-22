Ilastik Interface
=================


Ilastik Installation
--------------------

- install ilastik version 0.5 from

https://github.com/ilastik/ilastik-0.5/archive/v0.5.12.tar.gz

overwrite files from
./Interface/Ilastik/Ilastik
to
PATH_TO_ILASTIK_0.5/modules/classification/core/classifiers/


Vigra Installation
------------------

- install vigra / vigranumpy


under linux:

download vigra from:

http://ukoethe.github.io/vigra/doc-release/vigra/Installation.html

cmake -DLIB_SUFFIX=64
make 
make install

check ldconfig / /etc/ld.so.conf to make sure /usr/local/lib64 is in the library path


h5py Installation
------------------

* install h5py

  under linux:
  download repository
  
* matlab uses its own outdata version lighdf5.so.6
  to make h5py in matlab 2014b work:

* locate libhdf5.so.XXX typically in /usr/lib64

* relink 
  /usr/local/MATLAB/R2014b/sys/os/glnxa64/libhdf5.so.6 
  to the libhdf5.so.XXX, e.g.
  /usr/lib64/libhdf5.so.8

  and
  libhdf5_hl.so.6 
  to 
  /usr/lib64/libhdf5_hl.so.8



Ilastik under Matlab
--------------------

    initialize
    ilinitialize

then use

    ilclassify

see

    ./Interface/Ilastik/Test/testIlastik.m

for exmaples













