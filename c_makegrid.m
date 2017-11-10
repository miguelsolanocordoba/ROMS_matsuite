%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Build a ROMS grid file
%
%  Further Information:  
%  http://www.brest.ird.fr/Roms_tools/
%  
%  This file is part of ROMSTOOLS
%
%  ROMSTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  ROMSTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2002-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%  Contributions of P. Marchesiello (IRD) and X. Capet (UCLA)
%
%  Updated    Aug-2006 by Pierrick Penven
%  Updated    24-Oct-2006 by Pierrick Penven (mask correction)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; path(pathdef);
%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%
rpp_param
%
%%%%%%%%%%%%%%%%%%%%%% PATHS needed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath utility/mask
addpath utility/preprocessing_tools
addpath utility/mexcdf/netcdf_toolbox/netcdf
addpath utility/mexcdf/netcdf_toolbox/netcdf/ncsource
addpath utility/mexcdf/netcdf_toolbox/netcdf/nctype
addpath utility/mexcdf/netcdf_toolbox/netcdf/ncutility
addpath utility/mexcdf/mexnc
addpath utility/mexcdf/snctools
addpath utility/m_map
%%%%%%%%%%%%%%%%%%% END USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
warning off
%
% Title
%
disp(' ')
disp([' Making the grid: ',grdname])
disp(' ')
disp([' Title: ',ROMS_title])
disp(' ')
disp([' Resolution: 1/',num2str(1/dl),' deg'])
%
% Get the Longitude
%
lonr=(lonmin:dl:lonmax);
%
% Get the latitude for an isotropic grid
%
i=1;
latr(i)=latmin;
while latr(i)<=latmax
  i=i+1;
  latr(i)=latr(i-1)+dl*cos(latr(i-1)*pi/180);
end
[Lonr,Latr]=meshgrid(lonr,latr);
[Lonu,Lonv,Lonp]=rho2uvp(Lonr); 
[Latu,Latv,Latp]=rho2uvp(Latr);
%
% Create the grid file
%
disp(' ')
disp(' Create the grid file...')
[M,L]=size(Latp);
disp([' LLm = ',num2str(L-1)])
disp([' MMm = ',num2str(M-1)])
create_grid(L,M,grdname,ROMS_title)
%
% Fill the grid file
%
disp(' ')
disp(' Fill the grid file...')
nc=netcdf(grdname,'write');
nc{'lat_u'}(:)=Latu;
nc{'lon_u'}(:)=Lonu;
nc{'lat_v'}(:)=Latv;
nc{'lon_v'}(:)=Lonv;
nc{'lat_rho'}(:)=Latr;
nc{'lon_rho'}(:)=Lonr;
nc{'lat_psi'}(:)=Latp;
nc{'lon_psi'}(:)=Lonp;
close(nc)
%
%  Compute the metrics
%
disp(' ')
disp(' Compute the metrics...')
[pm,pn,dndx,dmde]=get_metrics(grdname);
xr=0.*pm;
yr=xr;
for i=1:L
  xr(:,i+1)=xr(:,i)+2./(pm(:,i+1)+pm(:,i));
end
for j=1:M
  yr(j+1,:)=yr(j,:)+2./(pn(j+1,:)+pn(j,:));
end
[xu,xv,xp]=rho2uvp(xr);
[yu,yv,yp]=rho2uvp(yr);
dx=1./pm;
dy=1./pn;
dxmax=max(max(dx/1000));
dxmin=min(min(dx/1000));
dymax=max(max(dy/1000));
dymin=min(min(dy/1000));
disp(' ')
disp([' Min dx=',num2str(dxmin),' km - Max dx=',num2str(dxmax),' km'])
disp([' Min dy=',num2str(dymin),' km - Max dy=',num2str(dymax),' km'])
%
%  Angle between XI-axis and the direction
%  to the EAST at RHO-points [radians].
%
angle=get_angle(Latu,Lonu);
%
%  Coriolis parameter
%
f=4*pi*sin(pi*Latr/180)/(24*3600);
%
% Fill the grid file
%
disp(' ')
disp(' Fill the grid file...')
nc=netcdf(grdname,'write');
nc{'pm'}(:)=pm;
nc{'pn'}(:)=pn;
nc{'dndx'}(:)=dndx;
nc{'dmde'}(:)=dmde;
nc{'x_u'}(:)=xu;
nc{'y_u'}(:)=yu;
nc{'x_v'}(:)=xv;
nc{'y_v'}(:)=yv;
nc{'x_rho'}(:)=xr;
nc{'y_rho'}(:)=yr;
nc{'x_psi'}(:)=xp;
nc{'y_psi'}(:)=yp;
nc{'angle'}(:)=angle;
nc{'f'}(:)=f;
nc{'spherical'}(:)='T';
close(nc);
%
%
%  Add topography from topofile
%
disp(' ')
disp(' Add topography...')
switch topofile
  case '~/data/etopo1.nc'
     h_roms=add_topo_ETOPO(grdname,topofile);
  case '~/data/SRTM30_PLUS_w100n40.nc'
     h_roms=add_topo_SRTM(grdname,topofile);
  case '~/data/prvi_1s_mhw_flip_roms.nc'
     h_roms=add_topo_DEM(grdname,topofile);
  otherwise
     h_roms=add_topo(grdname,topofile);
end

%% Interpolate raw bathymetry (for reference only) 
%% Data-set 
x=ncread(topofile,'x');
y=ncread(topofile,'y');
z=-(ncread(topofile,'z'));

[xx,yy]=meshgrid(x',y'); 
hraw=interp2(xx,yy,z',Lonr,Latr);

%save(['h_raw_' ROMS_config '.mat'],'hraw'); 

%
% Compute the mask
%
maskr=h_roms>0;
maskr=process_mask(maskr);
[masku,maskv,maskp]=uvp_mask(maskr);
%
%  Write it down
%
nc=netcdf(grdname,'write');
nc{'h'}(:)=h_roms;
nc{'mask_u'}(:)=masku;
nc{'mask_v'}(:)=maskv;
nc{'mask_psi'}(:)=maskp;
nc{'mask_rho'}(:)=maskr;
close(nc);
%
% Create the coastline
%

%if ~isempty(coastfileplot)
%  make_coast(grdname,coastfileplot);
%end

%
r=input('Do you want to use editmask ? y,[n]','s');
if strcmp(r,'y')
  disp(' Editmask:')
  disp(' Edit manually the land mask.')
  disp(' Press enter when finished.')
  disp(' ')
  if ~isempty(coastfileplot)
    editmask(grdname,coastfilemask)
  else
    editmask(grdname)
  end
  r=input(' Finished with edit mask ? [press enter when finished]','s');
end
%
close all
%
%  Smooth the topography
%
nc=netcdf(grdname,'write');
h=nc{'h'}(:);
maskr=nc{'mask_rho'}(:);
%
h=smoothgrid(h,maskr,hmin,hmax_coast,hmax,...
             rtarget,n_filter_deep_topo,n_filter_final);
%
%  Write it down
%
disp(' ')
disp(' Write it down...')
nc{'h'}(:)=h;
close(nc);
%
% make a plot
%
ld_library_path=getenv('DYLD_LIBRARY_PATH');
ld_library_path=['/opt/local/lib:/usr/local/lib' ld_library_path];
setenv('DYLD_LIBRARY_PATH',ld_library_path);

load ~/data/caribe_coast.mat

figure('visible','off')
themask=ones(size(maskr));
themask(maskr==0)=NaN; 
pcolor(Lonr,Latr,(-h).*themask);shading interp
hold on
plot(xcoast,ycoast,'k'); 
xlabel('longitude')
ylabel('latitude')
hold on
%[C1,h1]=contour(Lonr,Latr,h,[hmin 100 200 500 1000 2000 4000],'k');
%clabel(C1,h1,'LabelSpacing',1000,'Rotation',0,'Color','r')
xfac=cos(mean(mean(18.5))*pi/180); daspect([1 xfac 1]);
colorbar
caxis([-5000 0])  
print('-dpng','-r325',['bathy_' ROMS_config '.png'])
%0!convert -density 600 bathy.ps bathy.png

figure('visible','off')
themask=ones(size(maskr));
themask(maskr==0)=NaN;
pcolor(Lonr,Latr,(-h+hraw).*themask);shading interp
hold on
plot(xcoast,ycoast,'k');
xlabel('longitude')
ylabel('latitude')
hold on
%[C1,h1]=contour(Lonr,Latr,h,[hmin 100 200 500 1000 2000 4000],'k');
%clabel(C1,h1,'LabelSpacing',1000,'Rotation',0,'Color','r')
xfac=cos(mean(mean(18.5))*pi/180); daspect([1 xfac 1]);
colorbar
%caxis([-5000 0])
print('-dpng','-r325',['bathy_bias_' ROMS_config '.png'])

warning on
%!rm *.mat
%
% End
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
