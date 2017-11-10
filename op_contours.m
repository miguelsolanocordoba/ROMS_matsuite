clear all; close all; clc
%% OP_CONTOURS plot PROROMS contours
% 
% 

%------------------------------------ INPUT --------------------------------------%
wopt=input('\nPlot from thredds?y[n]: ','s')
if wopt=='y'
   valdate=input('\nEnter PROROMS simulation date ("yyyymmdd"): ','s')
else 
   wpath=input('\nWrite the path of history files to plot: ','s')
end

wcont=input('\nWrite contours to plot (e.g. velocity, temperature, salinity): ','s') 
wdays=input('\nEnter days to plot [e.g. 234]: ','s') 
wout=input('\nEnter output directory location: ','s')

%----------------------------------------------------------------------------------%

% Sort out days to plot and make file name strings 
numdays=size(wdays,2); 
for i=1:numdays
    days(i)=str2num(wdays(i)); 
end

count=0; 
switch wopt
   case 'y'
      for i=(1+24*(days(1)-1)):(numdays*24+1)
          count=count+1; 
          num=num2str(sprintf('%.4d',i)); 
          rfile(count,:)=['http://dm2.caricoos.org/thredds/dodsC/content/roms/' valdate '/r1/ocean_his_' num '.nc'];
          rtime(count)=addtodate(romstime(rfile(count,:)),-4,'hour'); 
      end
   otherwise
      for i=(1+24*(days(1)-1)):(numdays*24+1)
          count=count+1; 
          num=num2str(sprintf('%.4d',i)); 
          rfile(count,:)=[wpath '/ocean_his_' num '.nc']; 
          rtime(count)=addtodate(romstime(rfile(count,:)),-4,'hour'); 
      end
end

% Get grid info 
%datainfo=ncinfo(rfile(1,:));
temp=ncread(rfile(1,:),'u_eastward'); 
[xi_rho,eta_rho,N]=size(temp); 
lon_rho=ncread(rfile(1,:),'lon_rho');
lat_rho=ncread(rfile(1,:),'lat_rho');

% Load coast data
load caribe_coast

% 
fprintf('\nMaking plots...\n') 

%% MAIN LOOP 
for i=1:size(rfile,1)
    fprintf('%d/%d: %s\n',i,size(rfile,1),rfile(i,:))

    switch wcont  
        case 'velocity' 
            u=ncread(rfile(i,:),'u_eastward',[1 1 N 1],[xi_rho eta_rho 1 1]);
            v=ncread(rfile(i,:),'v_northward',[1 1 N 1],[xi_rho eta_rho 1 1]);
            mag_roms=sqrt(u.^2.+v.^2.);

            p=8;
            fh=figure;
            set(fh, 'Visible', 'off');
            pcolor(lon_rho,lat_rho,mag_roms);shading interp
            hold on
            quiver(lon_rho(1:p:end,1:p:end),lat_rho(1:p:end,1:p:end),...
            u(1:p:end,1:p:end),v(1:p:end,1:p:end),'k')
            hold on
            plot(xcoast,ycoast,'k')
            caxis([0 1]);
            hcol=colorbar;
            set(get(hcol,'title'),'String', '[m/s]');
            title(['ROMS ' datestr(rtime(i),0) ' local time']);
            xlabel('longitude');
            ylabel('latitude');
            xfac=cos(mean(mean(18.5))*pi/180); daspect([1 xfac 1]);
            print('-dpsc','-painter','-r325',['prvi_spd_roms_' num2str(sprintf('%.4d',i)) '.ps']);
            system(['convert -density  500 prvi_spd_roms_' num2str(sprintf('%.4d',i)) ...
                      '.ps prvi_spd_roms_' num2str(sprintf('%.4d',i)) '.png']);
     

        case 'temperature' 
            u=ncread(rfile(i,:),'u_eastward',[1 1 N 1],[xi_rho eta_rho 1 1]);
            v=ncread(rfile(i,:),'v_northward',[1 1 N 1],[xi_rho eta_rho 1 1]);
            t=ncread(rfile(i,:),'temp',[1 1 N 1],[xi_rho eta_rho 1 1]);
            mag_roms=sqrt(u.^2.+v.^2.);

            p=8;
            fh=figure;
            set(fh, 'Visible', 'off');
            pcolor(lon_rho,lat_rho,t);shading interp
            hold on
            quiver(lon_rho(1:p:end,1:p:end),lat_rho(1:p:end,1:p:end),...
            u(1:p:end,1:p:end),v(1:p:end,1:p:end),'k')
            hold on
            plot(xcoast,ycoast,'k')
            caxis([26 29]);
            hcol=colorbar;
            set(get(hcol,'title'),'String', '[m/s]');
            title(['ROMS ' datestr(rtime(i),0) ' local time']);
            xlabel('longitude');
            ylabel('latitude');
            xfac=cos(mean(mean(18.5))*pi/180); daspect([1 xfac 1]);
            print('-dpsc','-painter','-r325',['prvi_temp_roms_' num2str(sprintf('%.4d',i)) '.ps']);
            system(['convert -density  500 prvi_temp_roms_' num2str(sprintf('%.4d',i)) ...
                      '.ps prvi_temp_roms_' num2str(sprintf('%.4d',i)) '.png']);

        case 'salinity' 
            u=ncread(rfile(i,:),'u_eastward',[1 1 N 1],[xi_rho eta_rho 1 1]);
            v=ncread(rfile(i,:),'v_northward',[1 1 N 1],[xi_rho eta_rho 1 1]);
            s=ncread(rfile(i,:),'salt',[1 1 N 1],[xi_rho eta_rho 1 1]);
            mag_roms=sqrt(u.^2.+v.^2.);

            p=8;
            fh=figure;
            set(fh, 'Visible', 'off');
            pcolor(lon_rho,lat_rho,s);shading interp
            hold on
            quiver(lon_rho(1:p:end,1:p:end),lat_rho(1:p:end,1:p:end),...
            u(1:p:end,1:p:end),v(1:p:end,1:p:end),'k')
            hold on
            plot(xcoast,ycoast,'k')
            caxis([34.5 36.5]);
            hcol=colorbar;
            set(get(hcol,'title'),'String', '[m/s]');
            title(['ROMS ' datestr(rtime(i),0) ' local time']);
            xlabel('longitude');
            ylabel('latitude');
            xfac=cos(mean(mean(18.5))*pi/180); daspect([1 xfac 1]);
            print('-dpsc','-painter','-r325',['prvi_salt_roms_' num2str(sprintf('%.4d',i)) '.ps']);
            system(['convert -density  500 prvi_salt_roms_' num2str(sprintf('%.4d',i)) ...
                      '.ps prvi_salt_roms_' num2str(sprintf('%.4d',i)) '.png']);

    end
end
    
system(['rm *.ps']);
eval(['!mkdir ' wout]);
system(['mv *.png ' wout '/']);
