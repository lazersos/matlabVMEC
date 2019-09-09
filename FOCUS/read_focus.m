function data = read_focus( file )
%READ_FOCUS(file) Reads focus output files into a focus data structure
%   READ_FOCUS(file)  Reads the HDF5 file produced by the Focus code into a
%   data structure returned by the function.
%
%   Usage:
%   focus_data=read_focus('test.fo.h5');
%
%   See also read_hdf5.
%
%   Written by:     S. Lazerson (lazerson@pppl.gov)
%   Version:        1.0
%   Date:           02/08/17


data=read_hdf5(file);
data.file_ext=file(1:end-6);
data.rsurf=sqrt(data.xsurf.^2+data.ysurf.^2);
data.psurf=atan2(data.ysurf,data.xsurf);
data.psurf(data.psurf<0) = data.psurf(data.psurf<0)+2*pi;
data.datatype='focus_coil';

end

