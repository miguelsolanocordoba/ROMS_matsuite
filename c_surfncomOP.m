clear all; close all;path(pathdef);
%
disp('Starting surface ocnditions')
% Roms Pre Processing Parameters
%
rpp_param_new
%
% Extracting from AmSeas NCOM
ND=load('ncomdate.out','-ascii');
ncomdate=num2str(ND(1));
ncomdatefcstr=num2str(ND(2));
ncomdatefcnum=datenum(ncomdatefcstr,'yyyymmdd');

InterpMethod='regular'; 

counter=1;
btcounter=1;
for i=1:numfields

    try
      if i>8  %% Days 2-5
         ncomdatenum=addtodate(ncomdatefcnum,3*(i-1),'hour');
         num=sprintf('%.3d',3*(i-9));

         data=['http://ecowatch.ncddc.noaa.gov/thredds/dodsC/amseas/'...
          'ncom_relo_amseas_u_' ncomdate '00_t' num '.nc']

         fprintf('\nInterpolating surface conditions from AmSeas NCOM: %s',...
                  datestr(ncomdatenum,31));

       switch InterpMethod
         case 'regular' 
           [~,~,sst(:,:,counter),...
            shflux(:,:,counter),swflux(:,:,counter),swrad(:,:,counter)]=...
            extncomsurfOP(data,grdname,ncomdatenum);
         case 'scatter'
           [~,~,sst(:,:,counter),...
            shflux(:,:,counter),swflux(:,:,counter),swrad(:,:,counter)]=...
            extncomwindscatt(data,Gout);
        end

        sms_time(counter)=3*(btcounter-1)*60*60;


      else
         ncomdatenum=addtodate(ncomdatefcnum,3*(i-1),'hour');
         num=sprintf('%.3d',3*(i-1));

         data=['http://ecowatch.ncddc.noaa.gov/thredds/dodsC/amseas/'...
          'ncom_relo_amseas_u_' ncomdatefcstr '00_t' num '.nc']

         fprintf('\nInterpolating surface conditions from AmSeas NCOM: %s',...
                  datestr(ncomdatenum,31));

       switch InterpMethod
         case 'regular'
           [~,~,sst(:,:,counter),...
            shflux(:,:,counter),swflux(:,:,counter),swrad(:,:,counter)]=...
            extncomsurfOP(data,grdname,ncomdatenum);
         case 'scatter'
           [~,~,sst(:,:,counter),...
            shflux(:,:,counter),swflux(:,:,counter),swrad(:,:,counter)]=...
            extncomwindscatt(data,Gout);
        end

        sms_time(counter)=3*(btcounter-1)*60*60;

      end
        counter=counter+1;
        btcounter=btcounter+1;

    catch err
        disp('brinco 3 horas')
        disp(err)
        btcounter=btcounter+1;
        continue
    end
end

%%
grdinfo=ncinfo(grdname);

for i=1:length([grdinfo.Dimensions.Length])
    vv=genvarname(grdinfo.Dimensions(i).Name);
    eval([vv ' =grdinfo.Dimensions(i).Length;']);
end


system(['rm ' blkname]);
frcfile=blkname;

nccreate(frcfile,'SST','Dimensions',{'xi_rho',xi_rho,'eta_rho',eta_rho,...
    'sms_time',Inf},'Format','64bit')
ncwrite(frcfile,'SST',sst)
ncwriteatt(frcfile,'SST','long_name','sea surface temperature climatology')
ncwriteatt(frcfile,'SST','units','Celsius')
ncwriteatt(frcfile,'SST','time','sms_time')

nccreate(frcfile,'swrad','Dimensions',{'xi_rho',xi_rho,'eta_rho',eta_rho,...
    'sms_time',Inf},'Format','64bit')
ncwrite(frcfile,'swrad',swrad);
ncwriteatt(frcfile,'swrad','long_name','solar shortwave radiation flux');
ncwriteatt(frcfile,'swrad','units','Celsius meter second-1');
ncwriteatt(frcfile,'swrad','positive_value','downward flux, heating');
ncwriteatt(frcfile,'swrad','negative_value','upward flux, cooling');
ncwriteatt(frcfile,'swrad','time','sms_time');

nccreate(frcfile,'swflux','Dimensions',{'xi_rho',xi_rho,'eta_rho',eta_rho,...
    'sms_time',Inf},'Format','64bit');
ncwrite(frcfile,'swflux',swflux);
ncwriteatt(frcfile,'swflux','long_name','surface net freswater flux, (E-P)');
ncwriteatt(frcfile,'swflux','units','meter second-1');
ncwriteatt(frcfile,'swflux','time','sms_time');

nccreate(frcfile,'shflux','Dimensions',{'xi_rho',xi_rho,'eta_rho',eta_rho,...
    'sms_time',Inf},'Format','64bit');
ncwrite(frcfile,'shflux',shflux);
ncwriteatt(frcfile,'shflux','long_name','surface net heat flux')
ncwriteatt(frcfile,'shflux','units','Celsius meter second-1');
ncwriteatt(frcfile,'shflux','time','sms_time');

nccreate(frcfile,'sms_time','Dimensions',{'sms_time' Inf},'Format','64bit')
ncwrite(frcfile,'sms_time',sms_time)
ncwriteatt(frcfile,'sms_time','long_name','sea surface observation time')
ncwriteatt(frcfile,'sms_time','units','seconds')

ncwriteatt(frcfile,'/','Title','Wind Forcing netcdf file');
ncwriteatt(frcfile,'/','Base_Date',ncomdatefcstr);
ncwriteatt(frcfile,'/','creation_date',datestr(now));
