load alldates.txt  %yyyymmdd numbers
pth=28
ndates=length(alldates);
slcdir  = 'slcs/';

%polygon stuff for later - move to other script?
fid=fopen('latlon.txt','r');
a=fgetl(fid);
fclose(fid);

b=strsplit(a(2:end-1),','); %get rid of brackets
minlat=str2num(b{1});
maxlat=str2num(b{2});
minlon=str2num(b{3});
maxlon=str2num(b{4});

jnk=[minlon minlat maxlon minlat maxlon maxlat minlon maxlat minlon minlat];
polygon='polygon=';
for i=1:9
    polygon=[polygon num2str(jnk(i)) ','];
end
polygon=[polygon num2str(jnk(end))];

rootadr='https://api.daac.asf.alaska.edu/services/search/param?';

for i=1:ndates
    dn(i)=datenum(num2str(alldates(i)),'yyyymmdd');
end

for i=2:ndates
    
    %check if slc already coregistered
    flatslc    = [slcdir '/' num2str(alldates(i)) '.slc.full'];
    if(exist(flatslc))
        disp([num2str(alldates(i)) ' already made and flattened']);
    else
        disp(['working on ' num2str(alldates(i))]);
        %search for raw slc data files
        
        apicall=[rootadr polygon '&platform=Sentinel-1a,Sentinel-1b&processingLevel=SLC&start=' datestr(dn(i)) '&end=' datestr(dn(i)+1) '&relativeOrbit=' num2str(pth) '&output=metalink'];
        a=webread(apicall);
        b=regexp(a,'file name="(.{67}.zip)"','tokens');
        nfiles=length(b);
        disp(['found ' num2str(nfiles) ' files for date ' datestr(dn(i))]);
        
        %download raw slc data files
        
        for j=1:nfiles
            file=b{j}{1};
            if(exist(['data/' file]))
                disp([file ' already downloaded'])
            else
                command=['wget -c --http-user=rlohman --http-password=Rolohman1 https://datapool.asf.alaska.edu/SLC/SA/' file];
                system(command);
                movefile(file,'data'); %move to data dir
            end
        end
        
        
        %write date and int xml files and sequence of shell commands
        prep_datexml
        prep_intxml([1,2,3],[minlat maxlat minlon maxlon])
        %run
        system('./run_ints');
        %if slc file is complete, delete secondary raw data file (!!!!)
        if(exist(flatslc))
            for j=1:nfiles
                file=b{j}{1};
                delete(['data/' file]);
            end
            system(['rm -rf ' num2str(alldates(1)) '/int_' num2str(alldates(1)) '_' num2str(alldates(i))]);
        end
        
        
    end
end
