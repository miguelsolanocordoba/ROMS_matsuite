clear all; 
%
disp('Staring Boundary conditions')
% Roms Pre-Processing Parameters
%
rpp_param_new
%

%% Extract from AmSeas NCOM
ND=load('ncomdate.out','-ascii');
ncomdate=num2str(ND(1));
ncomdatefcstr=num2str(ND(2));
ncomdatefcnum=datenum(ncomdatefcstr,'yyyymmdd'); 

counter=1;
btcounter=1;
for i=1:numfields

    try
      if i>8  %% Days 2-5
         ncomdatenum=addtodate(ncomdatefcnum,3*(i-1),'hour');
         num=sprintf('%.3d',3*(i-9)); 

         data=['http://ecowatch.ncddc.noaa.gov/thredds/dodsC/amseas/'...
          'ncom_relo_amseas_u_' ncomdate '00_t' num '.nc'] 

         fprintf('\nInterpolating OBCs from AmSeas NCOM: %s',datestr(ncomdatenum,31)); 

         [uwest(:,:,counter),unorth(:,:,counter),ueast(:,:,counter),usouth(:,:,counter),...
         vwest(:,:,counter),vnorth(:,:,counter),veast(:,:,counter),vsouth(:,:,counter),...
         saltwest(:,:,counter),saltnorth(:,:,counter),salteast(:,:,counter),saltsouth(:,:,counter),...
         tempwest(:,:,counter),tempnorth(:,:,counter),tempeast(:,:,counter),tempsouth(:,:,counter),...
         surf_elwest(:,counter),surf_elnorth(:,counter),surf_eleast(:,counter),...
         surf_elsouth(:,counter)]=extncombry_newOP(data,Gout,ncomdatenum);
         
         bry_time(counter)=3*(btcounter-1)*60*60;

      else  %% Day 1
         ncomdatenum=addtodate(ncomdatefcnum,3*(i-1),'hour');
         num=sprintf('%.3d',3*(i-1));

         data=['http://ecowatch.ncddc.noaa.gov/thredds/dodsC/amseas/'...
          'ncom_relo_amseas_u_' ncomdatefcstr '00_t' num '.nc']

         fprintf('\nInterpolating OBCs from AmSeas NCOM: %s',datestr(ncomdatenum,31));

         [uwest(:,:,counter),unorth(:,:,counter),ueast(:,:,counter),usouth(:,:,counter),...
         vwest(:,:,counter),vnorth(:,:,counter),veast(:,:,counter),vsouth(:,:,counter),...
         saltwest(:,:,counter),saltnorth(:,:,counter),salteast(:,:,counter),saltsouth(:,:,counter),...
         tempwest(:,:,counter),tempnorth(:,:,counter),tempeast(:,:,counter),tempsouth(:,:,counter),...
         surf_elwest(:,counter),surf_elnorth(:,counter),surf_eleast(:,counter),...
         surf_elsouth(:,counter)]=extncombry_newOP(data,Gout,ncomdatenum);
      
         bry_time(counter)=3*(btcounter-1)*60*60;
       end
       
       counter=counter+1;
       btcounter=btcounter+1;
       
    catch err
       fprintf('\nError when interpolating field: %s',datestr(ncomdatenum,31));
       disp(err)
       btcounter=btcounter+1;
       continue
    end
end

%% Compute UBAR/VBAR
u.west=uwest;
u.east=ueast; 
u.north=unorth; 
u.south=usouth;
v.west=vwest;
v.east=veast; 
v.north=vnorth; 
v.south=vsouth;

[ubar,vbar]=uv_barotropic(u,v,Gout.Hz,[1 1 1 1]);

%% Create and writing Netcdf File

grdinfo=ncinfo(grdname);


for i=1:length([grdinfo.Dimensions.Length])
    vv=genvarname(grdinfo.Dimensions(i).Name);
    eval([vv ' =grdinfo.Dimensions(i).Length;']);
end

% importante verificar los siguientes parametros

% Input file dimensions
S.ncname=bryname;              % NetCDF file name
S.boundary=[1 1 1 1];          % Boundary switch to process 'west','east','south','north'
S.spherical=Gout.spherical;    % Spherical grid switch
S.Vtransform=Gout.Vtransform;  % Vertical transformation equation
S.Lm=Gout.Lm;                  % Number of interior RHO-points in X
S.Mm=Gout.Mm;                  % Number of interior RHO-points in Y
S.N=Gout.N;                    % Number of vertical levels
S.NT=2;                        % Number of active and passive tracer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c_boundaryfile(S,ncomdatefcnum);

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
ncwrite(S.ncname,'lon_rho_west',Gout.lon_rho_west);
ncwrite(S.ncname,'lat_rho_west',Gout.lat_rho_west);
ncwrite(S.ncname,'lon_rho_east',Gout.lon_rho_east);
ncwrite(S.ncname,'lat_rho_east',Gout.lat_rho_east);
ncwrite(S.ncname,'lon_rho_north',Gout.lon_rho_north);
ncwrite(S.ncname,'lat_rho_north',Gout.lat_rho_north);
ncwrite(S.ncname,'lon_rho_south',Gout.lon_rho_south);
ncwrite(S.ncname,'lat_rho_south',Gout.lat_rho_south);
ncwrite(S.ncname,'lon_u_west',Gout.lon_u_west);
ncwrite(S.ncname,'lat_u_west',Gout.lat_u_west);
ncwrite(S.ncname,'lon_u_east',Gout.lon_u_east);
ncwrite(S.ncname,'lat_u_east',Gout.lat_u_east);
ncwrite(S.ncname,'lon_u_north',Gout.lon_u_north);
ncwrite(S.ncname,'lat_u_north',Gout.lat_u_north);
ncwrite(S.ncname,'lon_u_south',Gout.lon_u_south);
ncwrite(S.ncname,'lat_u_south',Gout.lat_u_south);
ncwrite(S.ncname,'lon_v_west',Gout.lon_v_west);
ncwrite(S.ncname,'lat_v_west',Gout.lat_v_west);
ncwrite(S.ncname,'lon_v_east',Gout.lon_v_east);
ncwrite(S.ncname,'lat_v_east',Gout.lat_v_east);
ncwrite(S.ncname,'lon_v_north',Gout.lon_v_north);
ncwrite(S.ncname,'lat_v_north',Gout.lat_v_north);
ncwrite(S.ncname,'lon_v_south',Gout.lon_v_south);
ncwrite(S.ncname,'lat_v_south',Gout.lat_v_south);
ncwrite(S.ncname,'bry_time',bry_time);
ncwrite(S.ncname,'zeta_west',surf_elwest);
ncwrite(S.ncname,'zeta_east',surf_eleast);
ncwrite(S.ncname,'zeta_north',surf_elnorth);
ncwrite(S.ncname,'zeta_south',surf_elsouth);
ncwrite(S.ncname,'ubar_west',ubar.west);
ncwrite(S.ncname,'ubar_east',ubar.east);
ncwrite(S.ncname,'ubar_north',ubar.north);
ncwrite(S.ncname,'ubar_south',ubar.south);
ncwrite(S.ncname,'vbar_west',vbar.west);
ncwrite(S.ncname,'vbar_east',vbar.east);
ncwrite(S.ncname,'vbar_north',vbar.north);
ncwrite(S.ncname,'vbar_south',vbar.south);
ncwrite(S.ncname,'u_west',uwest);
ncwrite(S.ncname,'u_east',ueast);
ncwrite(S.ncname,'u_north',unorth);
ncwrite(S.ncname,'u_south',usouth);
ncwrite(S.ncname,'v_west',vwest);
ncwrite(S.ncname,'v_east',veast);
ncwrite(S.ncname,'v_north',vnorth);
ncwrite(S.ncname,'v_south',vsouth);
ncwrite(S.ncname,'temp_west',tempwest);
ncwrite(S.ncname,'temp_east',tempeast);
ncwrite(S.ncname,'temp_north',tempnorth);
ncwrite(S.ncname,'temp_south',tempsouth);
ncwrite(S.ncname,'salt_west',saltwest);
ncwrite(S.ncname,'salt_east',salteast);
ncwrite(S.ncname,'salt_north',saltnorth);
ncwrite(S.ncname,'salt_south',saltsouth);

ncwriteatt(S.ncname,'/','creation_time',datestr(now))
% ncwriteatt(S.ncname,'/','time_ref',datainfo.Attributes(14).Value);
