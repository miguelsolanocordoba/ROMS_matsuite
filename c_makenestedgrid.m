clear all;close all; %path(pathdef);
%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%
rpp_param4
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

addpath ~/Programs/PROROMS/matlab
addpath ~/Programs/PROROMS/matlab/landmask
addpath ~/Programs/PROROMS/matlab/utility
addpath ~/Programs/PROROMS/matlab/netcdf
addpath ~/Programs/PROROMS/matlab/bin
addpath ~/Programs/PROROMS/matlab/mex
addpath ~/Programs/PROROMS/matlab/grid
%%%%%%%%%%%%%%%%%%% END USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
warning off
%
%----------------------------- PART I: Parent Grid -----------------------------%
%
% Title
%
disp(' ')
disp([' Making the grid: ',grdname])
disp(' ')
disp([' Title: ',ROMS_title])
disp(' ')
disp([' Making Parent Grid ' ])
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
Lonu=Lonu';
Latu=Latu';
Lonv=Lonv';
Latv=Latv';
Lonp=Lonp';
Latp=Latp';
Lonr=Lonr';
Latr=Latr';
%
% Create the grid file
%
disp(' ')
disp(' Create the grid file...')
[Lp,Mp]=size(Lonr);
disp(' ') 
disp('Interior grid dimensions: ')
disp([' LLm = ',num2str(Lp-2)])
disp([' MMm = ',num2str(Mp-2)])

%create_grid(L,M,grdname,ROMS_title)
c_grid(Lp,Mp,grdname,true,true); 
%
% Fill the grid file
%
ncwrite(grdname,'lon_rho',Lonr)
ncwrite(grdname,'lat_rho',Latr)
ncwrite(grdname,'lon_u',Lonu)
ncwrite(grdname,'lat_u',Latu)
ncwrite(grdname,'lon_psi',Lonp)
ncwrite(grdname,'lat_psi',Latp)
ncwrite(grdname,'lon_v',Lonv)
ncwrite(grdname,'lat_v',Latv)
%
%  Compute the metrics
%
disp(' ')
disp('Compute the metrics...')
G.spherical=1;
G.lon_rho=Lonr;
G.lat_rho=Latr;
G.lon_u=Lonu;
G.lat_u=Latu;
G.lon_v=Lonv;
G.lat_v=Latv;
G.uniform=0;
[pm,pn,dndx,dmde]=grid_metrics(G,true);

[L,M]=size(Latp)

xr=0.*pm;
yr=xr;
for i=1:L
  xr(i+1,:)=xr(i,:)+2./(pm(i+1,:)+pm(i,:));
end
for j=1:M
  yr(:,j+1)=yr(:,j)+2./(pn(:,j+1)+pn(:,j));
end

% Interpolate metrics to uvp points
[xu,xv,xp]=rho2uvp(xr');
[yu,yv,yp]=rho2uvp(yr');

xu=xu';
xv=xv';
xp=xp';
yu=yu';
yv=yv';
yp=yp';

dx=1./pm;
dy=1./pn;
dxmax=max(max(dx/1000));
dxmin=min(min(dx/1000));
dymax=max(max(dy/1000));
dymin=min(min(dy/1000));
disp(' ')
disp('Grid statistics: ')
disp([' Min dx=',num2str(dxmin),' km - Max dx=',num2str(dxmax),' km'])
disp([' Min dy=',num2str(dymin),' km - Max dy=',num2str(dymax),' km'])
%
%  Angle between XI-axis and the direction
%  to the EAST at RHO-points [radians].
%
angle=get_angle(Latu',Lonu')';
%
%  Coriolis parameter
%
f=4*pi*sin(pi*Latr/180)/(24*3600);
%
xl=max(xp(:))-min(xp(:));
el=max(yp(:))-min(yp(:));
% Fill the grid file
%
disp(' ')
disp('Fill the grid file...')
ncwrite(grdname,'x_rho',xr)
ncwrite(grdname,'y_rho',yr)
ncwrite(grdname,'x_u',xu)
ncwrite(grdname,'y_u',yu)
ncwrite(grdname,'x_v',xv)
ncwrite(grdname,'y_v',yv)
ncwrite(grdname,'x_psi',xp)
ncwrite(grdname,'y_psi',yp)
ncwrite(grdname,'pm',pm)
ncwrite(grdname,'pn',pn)
ncwrite(grdname,'dndx',dndx)
ncwrite(grdname,'dmde',dmde)
ncwrite(grdname,'angle',angle)
ncwrite(grdname,'f',f)
ncwrite(grdname,'spherical',1)
ncwrite(grdname,'el',el)
ncwrite(grdname,'xl',xl)

%
%
%  Add topography from topofile
%
disp(' ')
disp('Add topography...')
switch topofile
  case '~/data/etopo1.nc'
     h_parent=add_topo_ETOPO(grdname,topofile);
     h_parent=h_parent'; 
  case '~/data/SRTM30_PLUS_w100n40.nc'
     h_parent=add_topo_SRTM(grdname,topofile);
     h_parent=h_parent'; 
  case '~/data/prvi_1s_mhw_flip_roms.nc'
     h_parent=add_topo_DEM(grdname,topofile);
     h_parent=h_parent'; 
  otherwise
     h_parent=add_topo(grdname,topofile); 
     h_parent=h_parent'; 
end
%h_parent=addtopo(grdname)';
%
% Check the bathymetry
%
pcolor(h_parent'); shading interp; colorbar
title('Parent grid bathymetry') 
xlabel('\xi')
ylabel('\eta')
pause(3) 
close all
%
% Compute the mask
%
disp('Computing mask...')
maskr=h_parent>0;
%maskr=process_mask(maskr);
[masku,maskv,maskp]=uvp_mask(maskr');
masku=masku';
maskv=maskv';
maskp=maskp';
maskr=double(maskr);
%
%  Write it down
%
ncwrite(grdname,'hraw',h_parent)
ncwrite(grdname,'mask_rho',maskr)
ncwrite(grdname,'mask_u',masku)
ncwrite(grdname,'mask_v',maskv)
ncwrite(grdname,'mask_psi',maskp)
%
% Edit the mask
%
addpath utility/mask
addpath utility/preprocessing_tools
addpath utility/mexcdf/netcdf_toolbox/netcdf
addpath utility/mexcdf/netcdf_toolbox/netcdf/ncsource
addpath utility/mexcdf/netcdf_toolbox/netcdf/nctype
addpath utility/mexcdf/netcdf_toolbox/netcdf/ncutility
addpath utility/mexcdf/mexnc
addpath utility/mexcdf/snctools
addpath utility/m_map

r=input('Do you want to edit the mask manually? y,[n]','s');
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
% Write the mask
ncwrite(grdname,'mask_rho',maskr)
ncwrite(grdname,'mask_u',masku)
ncwrite(grdname,'mask_v',maskv)
ncwrite(grdname,'mask_psi',maskp)
%
%  Smooth the topography
%
ncwrite(grdname,'h',h_parent); 
%
hPsmooth=smoothgrid(h_parent,maskr,hmin,hmax_coast,hmax,...
             rtarget,n_filter_deep_topo,n_filter_final);
%
%pcolor(hPsmooth); shading interp; colorbar
%  Write it down
%
disp(' ')
disp(' Write it down...')
ncwrite(grdname,'h',hPsmooth); 

%---------------------------- PART II: Grid Refinement --------------------------%
addpath ~/Programs/PROROMS/matlab
addpath ~/Programs/PROROMS/matlab/landmask
addpath ~/Programs/PROROMS/matlab/utility
addpath ~/Programs/PROROMS/matlab/netcdf
addpath ~/Programs/PROROMS/matlab/bin
addpath ~/Programs/PROROMS/matlab/mex
addpath ~/Programs/PROROMS/matlab/grid

disp(' ') 
disp(' Done with parent grid!!!')
disp(' ') 
disp([' Making nested grid (refinement): ', nstgrdname]) 
disp(' ') 
disp(' Locating coordinate indices form parent grid...')

%% Find closest coordinates in parent grid
[~,Imax]=min(abs(Lonp(:,1)-lonmaxR));
[~,Imin]=min(abs(Lonp(:,1)-lonminR));
[~,Jmax]=min(abs(Latp(1,:)-latmaxR));
[~,Jmin]=min(abs(Latp(1,:)-latminR));

%% Interpolate from parent (coarse) grid to child (fine) grid 
F = coarse2fine(grdname,nstgrdname,gf,Imin,Imax,Jmin,Jmax);

%% Change child grid topography (from Topo file)
disp(' ')
disp(' Add topography...')
switch topofile
  case '~/data/etopo1.nc'
     hR=add_topo_ETOPO(nstgrdname,topofile);
  case '~/data/SRTM30_PLUS_w100n40.nc'
     hR=add_topo_SRTM(nstgrdname,topofile);
    disp(' ENTRO AL CHILD')
  case '~/data/prvi_1s_mhw_flip_roms.nc'
     hR=add_topo_DEM(nstgrdname,topofile);
  otherwise
     hR=add_topo(nstgrdname,topofile);
end
%hR=addtopo(nstgrdname)';
%pcolor(hR);

% Compute maks 
disp(' ')
disp(' Computing mask...')
maskrR=hR>0;
[maskuR,maskvR,maskpR]=uvp_mask(maskrR);
maskuR=double(maskuR');
maskvR=double(maskvR');
maskpR=double(maskpR');
maskrR=double(maskrR');

%% Write new bathymetry and mask
%ncwrite(nstgrdname,'hraw',hR)
ncwrite(nstgrdname,'mask_rho',maskrR)
ncwrite(nstgrdname,'mask_u',maskuR)
ncwrite(nstgrdname,'mask_v',maskvR)
ncwrite(nstgrdname,'mask_psi',maskpR)


addpath utility/mask
addpath utility/preprocessing_tools
addpath utility/mexcdf/netcdf_toolbox/netcdf
addpath utility/mexcdf/netcdf_toolbox/netcdf/ncsource
addpath utility/mexcdf/netcdf_toolbox/netcdf/nctype
addpath utility/mexcdf/netcdf_toolbox/netcdf/ncutility
addpath utility/mexcdf/mexnc
addpath utility/mexcdf/snctools
addpath utility/m_map

r=input('Do you want to edit the mask manually? y,[n]','s');
if strcmp(r,'y')
  disp(' Editmask:')
  disp(' Edit manually the land mask.')
  disp(' Press enter when finished.')
  disp(' ')
  if ~isempty(coastfileplot)
    editmask(nstgrdname,coastfilemask)
  else
    editmask(nstgrdname)
  end
  r=input(' Finished with edit mask ? [press enter when finished]','s');
end
%
close all

%% Write down the new mask
disp(' ') 
disp(' Write down the new mask...') 
ncwrite(nstgrdname,'mask_rho',maskrR)
ncwrite(nstgrdname,'mask_u',maskuR)
ncwrite(nstgrdname,'mask_v',maskvR)
ncwrite(nstgrdname,'mask_psi',maskpR)

% SMOOTH BATHYMETRY %
disp(' ') 
disp(' Smoothing child grid bathymetry...')
hCsmooth=smoothgrid(hR',maskrR,hminR,hmax_coast,hmax,...
             rtargetR,n_filter_deep_topo,n_filter_final);

% Bathymetry nudging
%% Parent-child bathymetry nudging 
disp(' ') 
disp(' Nuding bathymetry from child grid to parent grid...') 

% Read child grid coordinates from new file 
lon_rhoR=ncread(nstgrdname,'lon_rho'); 
lat_rhoR=ncread(nstgrdname,'lat_rho'); 

[y1,y2]=meshgrid(lonr,latr);
hPinterp=interp2(y1,y2,hPsmooth',lon_rhoR,lat_rhoR);

% Load Nudging Coefficients
load NudgeCoef.mat

% Interpolate Nudging Coefficients to create alpha function 
xd1=linspace(lonr(1),lonr(end),numel(NudgeCoef(:,1)));
xd2=linspace(latr(1),latr(end),numel(NudgeCoef(1,:)));
xr1=linspace(lonr(1),lonr(end),size(lon_rhoR,1));
xr2=linspace(latr(1),latr(end),size(lat_rhoR,2));

[xD,yD]=meshgrid(xd1,xd2);
[xR,yR]=meshgrid(xr1,xr2);

alpha = interp2(xD,yD,5*NudgeCoef',xR,yR)';

% Nudgin to NCOM bathymetry (h_child = alpha*h_parent + (1-alpha)*h_fine)
hCnudg=alpha.*hPinterp+(1-alpha).*hCsmooth;
%pcolor(abs(hCsmooth-hCnudg)); shading interp; colorbar

%% Write it down
disp(' ') 
disp(' Write down new bathymetry...')
ncwrite(nstgrdname,'h',hCnudg);

% Write coordinates
addpath ~/Programs/PROROMS/matlab/grid
addpath ~/Programs/PROROMS/matlab/netcdf
addpath ~/Programs/PROROMS/matlab/utility
addpath ~/Programs/PROROMS/matlab/landmask

% Create contact points file 
disp(' ') 
disp(' Creating contact points file...') 
disp(' ')
[S,G]=contact({grdname,nstgrdname},ngcname);
%[S,G]=contact({grdname,nstgrdname},ngcname,true,false,false);

return
%--------------------------------- PART III: PLOTS -----------------------------%
ld_library_path=getenv('DYLD_LIBRARY_PATH');
ld_library_path=['/opt/local/lib:/usr/local/lib' ld_library_path];
setenv('DYLD_LIBRARY_PATH',ld_library_path);

figure('visible','off')

themask=ones(size(maskr));
themask(maskr==0)=NaN; 
themask2=ones(size(maskrR));
themask2(maskrR==0)=NaN; 
pcolor(Lonr,Latr,hPsmooth.*themask); shading interp
hold on
pcolor(lon_rhoR,lat_rhoR,hCnudg.*themask2); shading interp 
xlabel('longitude')
ylabel('latitude')
[C1,h1]=contour(Lonr,Latr,hPsmooth,[hmin 100 200 500 1000 2000 4000],'k');
clabel(C1,h1,'LabelSpacing',1000,'Rotation',0,'Color','r')
xfac=cos(mean(mean(18.5))*pi/180); daspect([1 xfac 1]);
colorbar
print('-dpng','-r325',['bathy_' ROMS_config '.png'])

warning on
%
% End
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
