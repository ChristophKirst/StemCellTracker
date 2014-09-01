Class Organization and Hierarchichy
===================================


(maybe incorporate tiling routines in the base classes???)
some experiments could need separate tiling for each frame
others might use the same tiling for each time point ?
similar for background 
-> have a routine the checks for this on the actuial level then goes
up if not present
-> easiest way




Base                       Tiled Image Data                 Patterned Image Data             


Experiment                 Experiment                       Experiment < Plate
^   ^                      ^   ^                            ^     ^
|   |                      |   |                            |     |
|   TimeSeries             |   TimeSeries                   |     TimeSeries(array) 
|   ^        ^             |   ^        ^                   |     ^         ^ 
|   |        |             |   |        |                   |     |         |
Frame        Trajectory    Frame   Trajectory               Colony(array) Trajectory  
^                          ^                                ^
|                          |                                |
|                          Tile                             Tile
|                          ^                                ^
|                          |                                |
Object < Cell              Object < Cell                    Object < Cell






Single Image Data          Tiled Image Data                 Patterned Image Data             


Experiment                 TiledExperiment                  TiledExperiment < Plate
^   ^                      ^   ^                            ^     ^
|   |                      |   |                            |     |
|   TimeSeries             |   TimeSeries                   |     TimeSeries(array) 
|   ^        ^             |   ^        ^                   |     ^         ^ 
|   |        |             |   |        |                   |     |         |
Frame        Trajectory    TiledFrame   Trajectory          Colony(array) Trajectory  
^                          ^                                ^
|                          |                                |
|                          Tile                             Tile
|                          ^                                ^
|                          |                                |
Object < Cell              Object < Cell                    Object < Cell
 



Experiment
----------

organizes information and routines related to the experiment
it is used to handle acces to the image data from various file formats, microscopes, etc
the data is stored in the results field


TimeSeries
----------

functionality for time series data, not necessary for static images
stores an array of frames and trajectories


Frame
-----

results and methods for a single image
stores segmented objects
base class for matching/tracking, it provides all routines needed by the tracking software


Trajectory
----------

stores results from the tracking and provides slicing routines etc.


Object
------

singel object representing the relevant data from the segmentation
necessary for the tracking, e.g. position and intensity information




TiledExperiment
---------------

derived form Experiment with additional functionality for tiled images



TiledFrame
----------

























TiledExperiment
-----
extends the Experiment class with additional functionality for tiling,
and storing positional relations between images





bla





need classes: ?
Alignment 
Surface


