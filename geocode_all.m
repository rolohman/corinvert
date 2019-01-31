pol='VV';
dem='/data/rlohman/Sentinel/Chile/dem/demLat_S26_S23_Lon_W071_W067.dem.wgs84';
bbox='-25.86 -23.00 -70.7 -68.1';

geodir=['geo_' pol];

if(~exist(geodir,'dir'))
    mkdir(geodir)
end

resdir=['results_dates_' pol '/'];

filelist=dir([resdir '*cor']);
nfiles=length(filelist)
for i=1:nfiles
    file=filelist(i).name;
    if(or(exist([geodir '/' file '.geo'],'file'),exist([resdir file '.geo'],'file')))
        disp([file ' already geocoded'])
    else
        geocoder(file,pol,resdir,bbox,dem,geodir);
    end
end

