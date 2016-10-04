% Copyright (c) 2012, Ilker R Capoglu 
% This script reads an Angora field-value file and plots the data.

fclose('all');
clear all;
    
fieldvaluefilename = '/my_path/FieldValueFile_Ey_0_0.hd5';

dt = hdf5_read(fieldvaluefilename,'time_step');
t0 = hdf5_read(fieldvaluefilename,'initial_time_value');

waveform = hdf5_read(fieldvaluefilename,'field_values');
length_time = length(waveform);
t = t0+dt*(0:length_time-1);

figure;
plot(t,waveform);
