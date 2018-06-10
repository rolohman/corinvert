function prep_intxml(swaths,boxgeom)
%box=[minlat maxlat minlon maxlon]
decide_ints
if(exist('dem','dir'))
    usedem=1;
    tmp=dir('dem/*dem.wgs84');
    demfile=tmp(1).name; %need more checks in here, plus option to put in own name
    
else
    usedem=0;
end
for i=1:length(ints)
    d1=dates(id1(i)).name;
    d2=dates(id2(i)).name;
    intname = ['int_' d1 '_' d2];
    datedir = [d1 '/'];
    intdir  = [datedir intname '/'];
    if(~exist(datedir,'dir'))
        mkdir(datedir);
    end
    if(~exist(intdir))
        mkdir(intdir)
    end
    
    file_xml=[intdir intname '.xml'];
    if(~exist(file_xml,'file'))
    fid = fopen(file_xml,'w');
    
    fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?>\n');
    fprintf(fid,'<topsApp>\n');
    fprintf(fid,'  <component name="topsinsar">\n');
    fprintf(fid,'    <property name="Sensor name">SENTINEL1</property>\n');
    fprintf(fid,'    <component name="master">\n');
    fprintf(fid,'      <catalog>../../%s_tops.xml</catalog>\n',d1);
    fprintf(fid,'    </component>\n');
    fprintf(fid,'    <component name="slave">\n');
    fprintf(fid,'      <catalog>../../%s_tops.xml</catalog>\n',d2);
    fprintf(fid,'    </component>\n');
    fprintf(fid,'    <property name="do unwrap">False</property> \n');
    fprintf(fid,'    <property name="unwrapper name">snaphu_mcf</property>\n');
    if(usedem)
        fprintf(fid,'    <property name="demFilename">%s</property>\n',demfile);
    end
    fprintf(fid,'    <property name="filter strength">0.1</property>\n');
    fprintf(fid,'    <property name="swaths">[');
    for j=1:length(swaths)
        fprintf(fid,'%d ',swaths(j));
    end
    fprintf(fid,']</property>\n');
    fprintf(fid,'     <property name="range looks">2</property>\n');
    fprintf(fid,'     <property name="azimuth looks">1</property>\n');
%    if(pth==166)
%        fprintf(fid,'<property name="region of interest">[34.610,34.736,-116.9,-114.102]</property>\n');
%    else
if(1)
        if(exist('boxgeom','var'))
            fprintf(fid,'<property name="region of interest">[');
            for i=1:3
                fprintf(fid,'%f,',boxgeom(i));
            end
            fprintf(fid,'%f]</property>\n',boxgeom(4));
        end
end
        %    end
    if(exist('boxgeom','var'))
        
        fprintf(fid,'    <property name="geocode bounding box">[');
        for i=1:3
            fprintf(fid,'%f,',boxgeom(i));
        end
        fprintf(fid,'%f]</property>\n',boxgeom(4));
    end
    fprintf(fid,'</component>\n');
    fprintf(fid,'</topsApp>\n');
    fclose(fid);
    end
end

%now make ints not already made
fid2 = fopen('run_ints','w');
home=pwd;
for i=1:ni
    d1=dates(id1(i)).name;
    d2=dates(id2(i)).name;
    intname = ['int_' d1 '_' d2];
    datedir = [d1 '/'];
    intdir  = [datedir intname '/'];
    if(~exist([intdir 'merged/topophase.cor.geo'],'file'))
        fprintf(fid2,['cd ' intdir '\n']);
        fprintf(fid2,['ln -s ' home '/dem/*dem*wgs84* .\n']);
        fprintf(fid2,['topsApp.py ' intname '.xml --steps --end=mergebursts\n']);
        fprintf(fid2,'rm -r %s %s coarse_coreg coarse_interferogram coarse_offsets ESD fine_offsets geom*\n',d1,d2);
        fprintf(fid2,'cd ../..\n');
    end
end
fclose(fid2);
