clear;close all;

%
% ROMS pre-processing parameters
%
rpp_param
%
% PATH
%
addpath utility/tides_otps/
addpath utility/mexcdf/mexnc/
addpath utility/mexcdf/snctools/
addpath utility/t_tide_v1.3beta/

t=datenum(2014,1,1);
tp=datenum(2015,1,1);

setenv('GFORTRAN_STDIN_UNIT', '5') 
setenv('GFORTRAN_STDOUT_UNIT', '6') 
setenv('GFORTRAN_STDERR_UNIT', '0')

ld_library_path=getenv('DYLD_LIBRARY_PATH');
ld_library_path=['/usr/local/lib:/opt/local/lib:/usr/local/lib/gcc/5:/usr/lib64/:/usr/lib:/lib64:/lib:' ld_library_path];
setenv('DYLD_LIBRARY_PATH',ld_library_path);

otps2frc_v4(grdname,t,tp,frcname,'OTPSnc/DATA/Model_AO_atlas');
