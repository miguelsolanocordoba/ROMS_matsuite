clear all;close all

addpath ~/data/roms_utility/timerefutils

% ROMS Pre-processng parameters

rpp_param
%bryname='~/PROROMS/OP_2.1/ROMS_FILES_20150206/roms_bryPR21.nc'
%frcname='~/PROROMS/OP_2.1/ROMS_FILES_20150206/roms_frcPR21.nc'
%grdname='~/PROROMS/OP_2.1/ROMS_FILES_20150206/roms_grdPR21.nc'
%theta_s=5; 
%theta_b=0.1; 
%hc=1; 
%N=32; 
%
%ntilei=4;
%ntilej=4;
% TIME_REF in Modified Julia days units

%time_ref='1968-05-23 00:00:00';

% Start date

%time_start=[ncomdate(1:4) '-' ncomdate(5:6) '-' ncomdate(7:8)];

% DSTART
%daydiff=abs(mjd(str2num(time_ref(1:4)),str2num(time_ref(6:7)),str2num(time_ref(9:10)))-...
%            mjd(str2num(time_start(1:4)),str2num(time_start(6:7)),str2num(time_start(9:10))));
%
%daydiffm=abs(datenum(time_ref)-datenum(time_start));
%
%dstart=daydiff-1;

% TIDE_START
time_ref_str=ncreadatt(bryname,'bry_time','units')
time_ref_num=datenum(time_ref_str(15:end))


tide_base_date=ncreadatt(frcname,'/','base_date');
tide_base_date_num=datenum(tide_base_date(12:end))
ty=str2num(tide_base_date(12:15));
tm=str2num(tide_base_date(17:18));
td=str2num(tide_base_date(20:21));

%tidedaydiff=abs(mjd(str2num(time_ref(1:4)),str2num(time_ref(6:7)),str2num(time_ref(9:10)))-...
%            mjd(ty,tm,td));
%tidedaydiffm=datenum(tide_base_date(12:end))-addtodate(datenum(time_start),50,'minute');
%tidedaydiffm=addtodate(datenum(tide_base_date(12:end)),-40,'minute')-datenum(time_start);

%%%%% change on apr 27, 2015
%tide_start=tide_base_date_num-time_ref_num

% FIX SHIFT 
tide_start=addtodate(tide_base_date_num,0,'minute')-time_ref_num

% check Lm and Mm

grdinfo=ncinfo(grdname);
Lm=grdinfo.Dimensions(5).Length-2;
Mm=grdinfo.Dimensions(6).Length-2;

%nstinfo=ncinfo(nstgrdname);
%LmN=nstinfo.Dimensions(1).Length-2;
%MmN=nstinfo.Dimensions(5).Length-2;

TIME_REF=[datestr(time_ref_num,'yyyymmdd'),'.00'];

fid = fopen('datestart.out','w');
fprintf(fid,'%11.11s\n',TIME_REF);
fprintf(fid,'%6.6f\n',tide_start);
fprintf(fid,'%6.6f\n',theta_s);
fprintf(fid,'%6.6f\n',theta_b);
fprintf(fid,'%6.6f\n',hc);
fprintf(fid,'%6.0f\n',N);
fprintf(fid,'%6.0f\n',Lm);
fprintf(fid,'%6.0f\n',Mm);
%fprintf(fid,'%6.0f\n',LmN);
%fprintf(fid,'%6.0f\n',MmN);
fprintf(fid,'%6.0f\n',ntilei);
fprintf(fid,'%6.0f\n',ntilej);
fclose(fid);
