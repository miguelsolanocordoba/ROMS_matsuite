%clear all;close all;path(pathdef);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                 %
% c_inincom.m: This script extracts initial conditions via OpenDAP from the       %
%                AmSeas NCOM server, interpolates values to the ROMS grid and     %
%                writes the initial conditions file in NetCDF format.             %
%                                                                                 %
% Variables extracted: u, v, salt, temp, surf_el, ubar, vbar                      %
%                                                                                 %
% Variable description:                                                           % 
%                                                                                 %
% u(xi_u,eta_u,s_rho,ocean_time) = eastward velocity component. [m/s]             %
%                                                                                 %
% v(xi_v,eta_v,s_rho,ocean_time) = northward velocity component. [m/s]            %
%                                                                                 %
% salt(xi_rho,eta_rho,s_rho,ocean_time) = salinity. [g/kg]                        %
%                                                                                 %
% temp(xi_rho,eta_rho,s_rho,ocean_time) = potential temprature. [C]               %
%                                                                                 %
% ubar(xi_u,eta_u,ocean_time) = eastward barotropic velocity component. [m/s]     %
%                                                                                 %
% vbar(xi_v,eta_v,ocean_time) = northward barotropic velocity component. [m/s]    %
%                                                                                 %
% zeta(xi_rho,eta_rho,ocean_time) = sea surface height. [m]                       %
%                                                                                 %
% Created: 01-May-2013 by E. Garcia and M. Solano                                 %
%                                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Roms Pre-Processing parameters
%
rpp_param_new
%
%-----------------------------Extract from AmSeas NCOM-----------------------------%
%% Read initial conditions (currently) from ecowatch server. 
% Start/read counters and dates
ND=load('ncomdate.out','-ascii'); 
ncomdate=num2str(ND(1)); 
ncomdatefcstr=num2str(ND(2)); 
ncomdatefcnum=datenum(ncomdatefcstr,'yyyymmdd'); 

% Main loop (Try OpenDAP ecowatch thredds server) 
try
    data=['http://ecowatch.ncddc.noaa.gov/thredds/dodsC/amseas/'...
        'ncom_relo_amseas_u_' ncomdatefcstr '00_t000.nc']

    tic
    [u,v,salt,temp,surf_el]=...
        extncomini_newOP(data,Gout,ncomdatefcnum);
    toc
    
    ini_time=0;

catch err
    disp('*** Error making initial conditions ***')
    fid=fopen('runstatus.txt','a'); 
    fprintf(fid,'\n%s\nEco-watch server is missing the 1st field.',...
                datestr(now,'mmm dd, yyyy'));
    fclose(fid); 
    disp(err)
    return
end

%% Compute UBAR/VBAR
[ubar,vbar]=uv_barotropic(u,v,Gout.Hz); 

%-----------------------------------------------------------------------------------%
% Input file dimensions
S.ncname=ininame;     %NetCDF file name
S.spherical=Gout.spherical;        %Spherical grid switch
S.Vtransform=Gout.Vtransform;       %Vertical transformation equation
S.Lm=Gout.Lm;        %Number of interior RHO-points in X
S.Mm=Gout.Mm;       %Number of interior RHO-points in Y
S.N=Gout.N;                %Number of vertical levels
S.NT=2;               %Number of active and passive tracer

% Create the file (from template)
c_initialfile(S,ncomdatefcnum);

% Write to NetCDF file
ncwrite(S.ncname,'spherical',Gout.spherical);
ncwrite(S.ncname,'Vtransform',Gout.Vtransform);
ncwrite(S.ncname,'Vstretching',Gout.Vstretching);
ncwrite(S.ncname,'theta_s',Gout.theta_s);
ncwrite(S.ncname,'theta_b',Gout.theta_b);
ncwrite(S.ncname,'Tcline',Gout.hc);
ncwrite(S.ncname,'hc',Gout.hc);
ncwrite(S.ncname,'s_rho',Gout.s_rho);
ncwrite(S.ncname,'s_w',Gout.s_w);
ncwrite(S.ncname,'Cs_r',Gout.Cs_r);
ncwrite(S.ncname,'Cs_w',Gout.Cs_w);
ncwrite(S.ncname,'h',Gout.h);
ncwrite(S.ncname,'lon_rho',Gout.lon_rho);
ncwrite(S.ncname,'lat_rho',Gout.lat_rho);
ncwrite(S.ncname,'lon_u',Gout.lon_u);
ncwrite(S.ncname,'lat_u',Gout.lat_u);
ncwrite(S.ncname,'lon_v',Gout.lon_v);
ncwrite(S.ncname,'lat_v',Gout.lat_v);
ncwrite(S.ncname,'ocean_time',ini_time);
ncwrite(S.ncname,'zeta',surf_el);
ncwrite(S.ncname,'ubar',ubar);
ncwrite(S.ncname,'vbar',vbar);
ncwrite(S.ncname,'u',u);
ncwrite(S.ncname,'v',v);
ncwrite(S.ncname,'temp',temp);
ncwrite(S.ncname,'salt',salt);

% Write creation time
ncwriteatt(S.ncname,'/','creation_time',datestr(now))

%------------------------------------ EOF --------------------------------------%
