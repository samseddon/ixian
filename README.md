# ixian
Integral X-ray Image ANalysis is a python based command line run tool designed
to take Dectris imaged X-Ray diffraction data, and convert it into reciprocal
space for analysis and publication level plotting. 


# Command line flags

"-dataset" will run code from start, allowing data folder and experiment
identifier from the start

"-oscan" will populate a 3D QSpace from a 2D dectris and omega scan

"-plot" will plot three cubic Q compressions 

"-h" will run the last scan number ran again meaning you don't have to 
re-enter the scan number at the command line
