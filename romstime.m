function [date_roms] = romstime(romsfile)
%% [date_roms] = romstime(romsfile) creates roms date vector. 
% The function romstime takes the file 'romsfile' which is the netcdf
% history output produced by ROMS and creates a Matlab date vector in
% datenum format. It uses attribute data from romsfile to create the date
% vector from the ocean_time variable. 

time_roms=ncread(romsfile,'ocean_time');
start_roms=ncread(romsfile,'dstart');
time_roms_units=ncreadatt(romsfile,'ocean_time','units');
time_roms_ref_str=time_roms_units(15:33);
time_roms_ref_num=datenum(time_roms_ref_str,0);
date_roms=zeros(1,length(time_roms));
for i=1:length(time_roms)
   date_roms(i)=addtodate(addtodate(time_roms_ref_num,start_roms,'day'),time_roms(i),'second');
end