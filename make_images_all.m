decide_ints

nx_geo    = 7560;
ny_geo    = 6228;
masterdir = [home '/' dates(1).name '/int_' dates(1).name '_' dates(2).name '/'];
TSdir     = [home '/results_TS/'];
invdir    = [home '/results_dates/'];
mergeddir = [masterdir 'merged/']; %master interferogram dir, location of topophase.cor for geocoding
togeodir  = [home '/geotiffs/togeo/'];
dir2      = [home '/geotiffs/geofiles/']; %directory for the geocoded r4 files
finaldir  = [home '/geotiffs/TS/']; %directory for images
coldir    = '/home/rlohman/matlab/ISCE_scripts/'; %directory for colormap txt files
nullval   = -9999; %null value for gdal/earthengine
coltype   = 'color'; %'color' or bw'


if(~exist([home '/geotiffs'],'dir'))
    mkdir([home '/geotiffs']);
end
if(~exist(togeodir,'dir'))
    mkdir(togeodir);
end
if(~exist(dir2,'dir'))
    mkdir(dir2);
end
if(~exist(finaldir,'dir'))
    mkdir(finaldir);
end
if(~exist([finaldir 'b1/'],'dir'))
    mkdir([finaldir 'b1/']);
end

%make and geocode maskfile for masking blank edges in geocoded ints
edgemaskfile=[home '/geotiffs/edgemask.r4'];
geomaskfile=[home '/geotiffs/geoedgemask.r4'];
if(~exist(geomaskfile,'file'))
    if(~exist(edgemaskfile,'file'))
        fid=fopen(edgemaskfile,'w');
        fwrite(fid,ones(newnx,newny),'real*4');
        fclose(fid);
    end
    command=['mag_phs2rmg ' edgemaskfile ' ' edgemaskfile ' ' mergeddir 'topophase.cor ' num2str(newnx)];
    system(command);
    chdir(masterdir);
    system(['topsApp.py smallproc.xml --steps --dostep=geocode']);
    chdir(home);
    
    command=['rmg2mag_phs ' mergeddir 'topophase.cor.geo tmp ' geomaskfile ' ' num2str(nx_geo)];
    system(command);
end

if(0)
    %rdate  = {'20150325','20150807','20160625'};
    rdate  = {'20150325','20150807','20160625','20170310','20170513','20170606'};
    dnr    = datenum(rdate,'yyyymmdd');
    rdate  = rdate(dnr<max(dn));
    dnr    = dnr(dnr<max(dn));
    
    for i=1:length(rdate)
        if(~exist([togeodir rdate{i} '.mag'],'file'));
            fidmag=fopen([TSdir rdate{i} '.mag10'],'r');
            fidmag10=fopen([TSdir rdate{i} '.mag'],'r');
            fidtim=fopen([TSdir rdate{i} '.time'],'r');
            fidmag0=fopen([togeodir rdate{i} '.mag'],'w');
            fidmag100=fopen([togeodir rdate{i} '.mag10'],'w');
            fidtim0=fopen([togeodir rdate{i} '.time'],'w');
            for j=1:newny
                mag=fread(fidmag,newnx,'real*4');
                mag10=fread(fidmag10,newnx,'real*4');
                tim=fread(fidtim,newnx,'real*4');
                id=find(and(tim>50,mag<0.1));
                mag(id)=0;
                mag10(id)=0;
                tim(id)=0;
                id=find(tim>149);
                mag(id)=0;
                mag10(id)=0;
                tim(id)=0;
                id=find(mag==0);
                mag(id)=0;
                mag10(id)=0;
                tim(id)=0;
                mag(mag>1)=1;
                mag10(mag10>1)=1;
                fwrite(fidmag0,mag,'real*4');
                fwrite(fidmag100,mag10,'real*4');
                fwrite(fidtim0,tim,'real*4');
            end
            fclose('all');
        end
    end
    if(~exist([togeodir 'resn'],'file'))
        system(['ln -s ' TSdir 'resn  ' togeodir]);
    end
end
files=[dir([invdir 'rel*cor']);dir([invdir 'perm*'])];
%files=[dir([invdir 'rel*filt.cor'])];
for i=1:length(files)
    if(~exist([togeodir files(i).name],'file'));
        system(['ln -s ' invdir files(i).name ' ' togeodir]);
    end
end

if(~exist([togeodir 'c0.cor'],'file'))
    system(['ln -s ' invdir 'c0.cor ' togeodir]);
end
if(~exist([togeodir 'rms.cor'],'file'))
    system(['ln -s ' invdir 'rms.cor ' togeodir]);
end

%files=[dir([togeodir 'resn']);dir([togeodir 'perm*']);dir([togeodir '2*']);dir([togeodir 'c*']);dir([togeodir 'rel*']);];
files=dir([togeodir]);
for i=1:length(files)
    files(i).name
    if(files(i).bytes==newnx*newny*4)
        if(~exist([finaldir files(i).name '.tif'],'file'))
            if(~exist([dir2 files(i).name],'file'))
                command=['mag_phs2rmg ' togeodir files(i).name ' '  togeodir files(i).name ' ' mergeddir 'topophase.cor ' num2str(newnx)];
                system(command);
                chdir(masterdir);
                system(['topsApp.py smallproc.xml --steps --dostep=geocode']);
                chdir(home);
                
                command=['rmg2mag_phs ' mergeddir 'topophase.cor.geo tmp ' dir2 files(i).name ' ' num2str(nx_geo)];
                system(command);
                
                fid  = fopen([dir2 files(i).name],'r');
                fidm = fopen(geomaskfile,'r');
                fido = fopen('tmp','w');
                for j=1:ny_geo
                    tmp         = fread(fid,nx_geo,'real*4');
                    msk         = fread(fidm,nx_geo,'real*4');
                    tmp(msk==0) = nullval;
                    tmp(isnan(tmp))=nullval;
                    fwrite(fido,tmp,'real*4');
                end
                fclose('all');
                movefile('tmp',[dir2 files(i).name]);
                
                
                fid=fopen([dir2 files(i).name '.vrt'],'w');
                fprintf(fid,'<VRTDataset rasterXSize="%d" rasterYSize="%d">\n',nx_geo,ny_geo);
                fprintf(fid,'<SRS>EPSG:4326</SRS>\n');
                %need to change below for other lat/lon geocoded areas
                fprintf(fid,'<GeoTransform>-70.6, 0.0002777777777777778, 0.0, -24.0, 0.0, -0.0002777777777777778</GeoTransform>\n');
                fprintf(fid,'<VRTRasterBand band="1" dataType="Float32" subClass="VRTRawRasterBand">\n');
                fprintf(fid,'    <SourceFilename relativeToVRT="1">%s</SourceFilename>\n',files(i).name);
                fprintf(fid,'    <ByteOrder>LSB</ByteOrder>\n');
                fprintf(fid,'    <ColorInterp>Palette</ColorInterp>\n');
                fprintf(fid,'    <ImageOffset>0</ImageOffset>\n');
                fprintf(fid,'    <PixelOffset>4</PixelOffset>\n');
                fprintf(fid,'    <LineOffset>%d</LineOffset>\n',nx_geo*4);
                fprintf(fid,'    <NoDataValue>%d</NoDataValue>\n',nullval);
                fprintf(fid,'</VRTRasterBand>\n');
                fprintf(fid,'</VRTDataset>\n');
                fclose(fid);
                
            end
            if(0)
            if(~exist([finaldir files(i).name '.tif'],'file'))
                if(regexp(files(i).name,'rel'))
                    command=['gdaldem color-relief  ' dir2 files(i).name '.vrt ' coldir coltype '_mag.txt ' finaldir files(i).name '.tif -alpha'];
                elseif(regexp(files(i).name,'perm'))
                    command=['gdaldem color-relief  ' dir2 files(i).name '.vrt ' coldir coltype '_mag.txt ' finaldir files(i).name '.tif -alpha'];
                elseif(regexp(files(i).name,'mag'))
                    command=['gdaldem color-relief  ' dir2 files(i).name '.vrt ' coldir coltype '_mag.txt ' finaldir files(i).name '.tif -alpha'];
                elseif(regexp(files(i).name,'time'))
                    command=['gdaldem color-relief  ' dir2 files(i).name '.vrt ' coldir coltype '_time.txt ' finaldir files(i).name '.tif -alpha'];
                elseif(regexp(files(i).name,'resn'))
                    command=['gdaldem color-relief  ' dir2 files(i).name '.vrt ' coldir coltype '_mag.txt ' finaldir files(i).name '.tif -alpha'];
                elseif(regexp(files(i).name,'c0.cor'))
                    command=['gdaldem color-relief  ' dir2 files(i).name '.vrt ' coldir coltype '_mag.txt ' finaldir files(i).name '.tif -alpha'];
                else
                    disp([files(i).name ' name not matched']);
                    return
                end
                disp(command)
                system(command)
            end
            end
        end

        outname=regexprep(files(i).name,'\.','_');
        if(~exist([finaldir 'b1/' outname '_b1.tif'],'file')) %these files are used by google earthengine
            command=['gdal_translate ' dir2 files(i).name '.vrt ' finaldir 'b1/' outname '_b1.tif'];
            disp(command);
            system(command);
        end
    else
        disp(['file ' files(i).name ' wrong size: ' files(i).bytes]);
    end
end


%old info:
% convert out.ppm -transparent "rgb(0,0,144)" test.png
% convert out.ppm -transparent "rgb(0,0,255)" test.png
% gdal_translate -of VRT -a_srs EPSG:4326 -a_ullr -70.6 -24 -68.5 -25.73    test.png newtest.vrt
% gdal_translate -of vrt -expand rgba newtest.vrt temp.vrt
% gdal2tiles.py -p geodetic temp.vrt -k
%         if(~exist([finaldir files(i).name],'dir'))
%             command=['gdal2tiles.py ' finaldir files(i).name '.tif -k ' finaldir files(i).name];
%             disp(command)
%             system(command);
%         end