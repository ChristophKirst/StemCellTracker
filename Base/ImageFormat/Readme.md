ImageFormat
===========

routines to reshape and permute data arrays into various formats


Formats
-------

   - base format: 'XYZCT'

   - base format axis inversions: 'xyzct'

   - cell formats: 'UVW'

   - cell formats axis inversion uvw


Any combitnation of the formats can be used e.g. 'yXCTz'


Special Formats
---------------

   - 'matlab' = 'yXCZT'

   - 'imagej' = 'yXZCT'


Examples
--------
   
   - 'XY'    = X x Y maTrix = 'XY' (2D grayscale)
   - 'XYC'   = X x Y x C matrix (2D multi channel )
   - 'XYZ'   = X x Y x Z matrix (3D)
   - 'XYCZ'  = X x Y x C x Z matrix (3D multi channel image)
   - 'XYZC'  = X x Y x Z x C matrix (3D multi channel image)
   - 'XYZT'  = X x Y x Z x T matrix (4D grayscale)
   - 'XYCZT' = X x Y x C x Z x T matrix (4D multi channel image)
   - 'XYZCT' = X x Y x Z x C x T matrix (4D multi channel image, time last)
   - 'XYZCT' = X x Y x Z x T x C matrix (4D multi channel image, channel lasy)









