function new_data = read_hdr(hdr_filename,filename)
% Quick script to import the spectral data subset for Frankenberg. Change
% The paths for the new_data_filename and new_data_hdr_filename. This also
% requires the read_envihdr.m function in the search path.

% File names for radiance cube and corresponding header 
%new_data_hdr_filename = 'Subset_for_Frankenberg.hdr';
%new_data_filename = 'Subset_for_Frankenberg';

hdr_data = read_envihdr(hdr_filename);
%[filepath,filename,ext] = fileparts(hdr_filename);
%filename = join([filepath, filename],"/");

new_data=multibandread(filename,[hdr_data.lines,hdr_data.samples,hdr_data.bands],...
    strcat(hdr_data.format,'=>float32'),hdr_data.header_offset,...
    hdr_data.interleave,hdr_data.machine);
end