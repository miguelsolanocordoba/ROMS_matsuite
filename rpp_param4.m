%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% RPP_PARAM4: Input parameter file for all PROROMS4 pre-processing files
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%---------------------------- PATHS AND FILES -----------------------------%
%
%  ROMS title names and directories
%
ROMS_title  = 'PROROMS40';
ROMS_config = 'PROROMS40';
%
%  ROMSTOOLS directory
%
ROMSTOOLS_dir = './utility/';
%
%  Run directory
%
RUN_dir=[pwd,'/'];
%
%  ROMS input netcdf files directory
%
ROMS_files_dir=[RUN_dir,'ROMS_FILES/'];
%
%  Global data directory (etopo, coads, datasets download from ftp, etc..)
%
DATADIR=ROMSTOOLS_dir; 
%
eval(['!mkdir ',ROMS_files_dir])
%
% ROMS file names (grid, forcing, bulk, climatology, initial)
%
% Parent Grid
hisname=[ROMS_files_dir,'ocean_his_0001.nc'];
grdname=[ROMS_files_dir,'roms_grid_parent_1_30.nc'];
ngcname=[ROMS_files_dir,'NGC_PROROMS.nc'];
bryname=[ROMS_files_dir,'BRY_P_PROROMS.nc'];
bryname_detied=[ROMS_files_dir,'FBRY_P_PROROMS.nc'];
ininame=[ROMS_files_dir,'INI_P_PROROMS.nc'];
blkname=[ROMS_files_dir,'BLK_P_PROROMS.nc'];
wndname=[ROMS_files_dir,'WND_P_PROROMS.nc'];
clmname=[ROMS_files_dir,'CLM_P_PROROMS.nc'];
frcname=[ROMS_files_dir,'FRC_P_PROROMS.nc'];
% 
% Child Grid
nstgrdname=[ROMS_files_dir,'roms_grid_child_CR3.nc'];
nstininame=[ROMS_files_dir,'INI_C_PROROMS.nc'];
nstblkname=[ROMS_files_dir,'BLK_C_PROROMS.nc'];
nstwndname=[ROMS_files_dir,'WND_C_PROROMS.nc'];
nstclmname=[ROMS_files_dir,'CLM_C_PROROMS.nc'];
nstfrcname=[ROMS_files_dir,'FRC_C_PROROMS.nc'];
%
%
%---------------------------- GRID PARAMETERS -----------------------------%
% 
% Grid Refinement (Nesting) 
gf=3;        % Refinement coefficient (Must be 3,5 or 7)
%
% Parent Grid dimensions:
%
%lonmin = -71.00;   % Minimum longitude [degree east]
%lonmax = -61.00;   % Maximum longitude [degree east]
%latmin =  15.00;   % Minimum latitude  [degree north]
%latmax =  20.60;   % Maximum latitude  [degree north]
lonmin = -71.822592082366618;   % Minimum longitude [degree east]
lonmax = -61.003688732598619;   % Maximum longitude [degree east]
latmin =  14.999000014414207;   % Minimum latitude  [degree north]
latmax =  20.698400023062729;   % Maximum latitude  [degree north]
%
%% Nested Grid dimensions 
lonminR = -67.9697;;   % Minimum longitude [degree east]
lonmaxR = -64.0367;   % Maximum longitude [degree east]
latminR =  17.2000;   % Minimum latitude  [degree north]
latmaxR =  19.1000;   % Maximum latitude  [degree north]
%
% Parent Grid resolution [degree]
%
dl = 1/30;
%
% Number of vertical Levels %
N = 32;
%
%  Vertical grid parameters 
%
theta_s = 3.;
theta_b = 0.4;
hc      = 100.;
%
% Minimum depth at the shore [m] (depends on the resolution,
% rule of thumb: dl=1, hmin=300, dl=1/4, hmin=150, ...)
% This affect the filtering since it works on grad(h)/h.
%
hmin = 10;
hminR =10;
%
% Maximum depth at the shore [m] (to prevent the generation
% of too big walls along the coast)
%
hmax_coast = 500;
%
% Maximum depth [m] (cut the topography to prevent
% extrapolations below WOA data)
%
hmax = 5000;
%
%  Topography netcdf file name (ETOPO 2 or any other netcdf file
%  in the same format)
%
topofile = ['~/data/SRTM30_PLUS_w100n40.nc'];
%
%
%---------------------------- SMOOTHING PARAMETERS -----------------------------%
%
% Slope parameter (r=grad(h)/h) maximum value for topography smoothing
%
rtarget = 0.25;
rtargetR = 0.25;
%
% Order of Shapiro filter (order) and maximum number of iterations (maxit)
%
order=4; 
maxit=1500; 
%
% Number of pass of a selective filter to reduce the isolated
% seamounts on the deep ocean.
%
n_filter_deep_topo=4;
%
% Number of pass of a single hanning filter at the end of the
% smooting procedure to ensure that there is no 2DX noise in the 
% topography.
%
n_filter_final=2;
%
%  GSHSS user defined coastline (see m_map) 
%  XXX_f.mat    Full resolution data
%  XXX_h.mat    High resolution data
%  XXX_i.mat    Intermediate resolution data
%  XXX_l.mat    Low resolution data
%  XXX_c.mat    Crude resolution data
%
coastfileplot = 'coastline_f.mat';
coastfilemask = 'coastline_f_mask.mat';
%%%%%%%%%%%%%%%%%%%%%555%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                   %
%  Date of initialization                           %
%                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncomdate='20160401';
numdays=30; 
numfields=numdays*8+1;
%[~,ncomdate_tmp]=system('date +%Y%m%d');
%ncomdate=ncomdate_tmp(1:8)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                   %
%  Domain decomposition                             %
%                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ntilei=8;
ntilej=4;
%%
InterpMethod='regular'; 
%
interp_method = 'linear';           % Interpolation method: 'linear' or 'cubic'
%
makeplot     = 0;                 % 1: create a few graphics after each preprocessing step
%
