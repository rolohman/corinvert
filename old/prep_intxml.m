function prep_intxml(swaths,boxgeom)
%box=[minlat maxlat minlon maxlon]
decide_ints
fid_run = fopen('run_ints','w');
slcdir  = 'slcs/';
if(~exist(slcdir,'dir'))
    mkdir(slcdir)
end

if(exist('dem','dir'))
    usedem=1;
    tmp=dir('dem/*dem.wgs84');
    demfile=tmp(1).name; %need more checks in here, plus option to put in own name
else
    usedem=0;
end

for i=1:nd-1   
    d1=dates(id1(i)).name;
    d2=dates(id2(i)).name;
    intname = ['int_' d1 '_' d2];
    datedir = [d1 '/'];
    intdir  = [datedir intname '/'];
    primary0   = [intdir 'merged/master.slc.full'];
    primary    = [slcdir 'master.slc.full'];
    secondary  = [intdir 'merged/slave.slc.full'];
    flatint    = [intdir 'merged/topophase.flat.full'];
    ramp       = [intdir 'merged/ramp.full'];
    flatslc    = [slcdir d2 '.slc.full'];
    if(i==1)
        intname0=intname;
    end
    if(exist(flatslc))
        disp([d2 ' already made and flattened']);
    else
        
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
            fprintf(fid,'      <catalog>../../datexml/%s_tops.xml</catalog>\n',d1);
            fprintf(fid,'    </component>\n');
            fprintf(fid,'    <component name="slave">\n');
            fprintf(fid,'      <catalog>../../datexml/%s_tops.xml</catalog>\n',d2);
            fprintf(fid,'    </component>\n');
            if(usedem)
                fprintf(fid,'    <property name="demFilename">%s</property>\n',demfile);
            end
            fprintf(fid,'    <property name="swaths">[');
            for j=1:length(swaths)
                fprintf(fid,'%d ',swaths(j));
            end
            fprintf(fid,']</property>\n');
            if(exist('boxgeom','var'))
                fprintf(fid,'<property name="region of interest">[');
                for i=1:3
                    fprintf(fid,'%f,',boxgeom(i));
                end
                fprintf(fid,'%f]</property>\n',boxgeom(4));
            end
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
    
        
        fprintf(fid_run,['cd ' intdir '\n']);
        fprintf(fid_run,['ln -s ' pwd '/dem/*dem*wgs84* .\n']);
        if(i==1)
            fprintf(fid_run,['topsApp.py ' intname '.xml --steps\n']);
        else
        
        fprintf(fid_run,['topsApp.py ' intname '.xml --steps --end=verifyDEM\n']);
        fprintf(fid_run,['ln -s ../' intname0 '/geom_master/ .\n']);
        fprintf(fid_run,['cp ../' intname0 '/PICKLE/topo* PICKLE\n']);
        fprintf(fid_run,['sed -i ''s/' dates(id2(1)).name '/' d2 '/g'' PICKLE/topo.xml\n']);
        fprintf(fid_run,['topsApp.py ' intname '.xml --steps --start=subsetoverlaps --end=mergebursts\n']);
%         fprintf(fid_run,['topsApp.py ' intname '.xml --steps --start=subsetoverlaps --end=rangecoreg\n']);
%         fprintf(fid_run,['topsApp.py ' intname '.xml --steps --dostep=fineoffsets\n']);
%         fprintf(fid_run,['topsApp.py ' intname '.xml --steps --start=burstifg --end=mergebursts\n']);
         end
        
        fprintf(fid_run,'cd ../..\n');
        
        if(i==1) %should be made with first interferogram
            command = ['looks.py -i ' primary0 ' -o ' primary ' -r 1 -a 1\n'];
            fprintf(fid_run,command);
        end
        
        command = ['looks.py -i ' secondary ' -o ' secondary ' -r 1 -a 1\n'];
        fprintf(fid_run,command);
        command=['looks.py -i ' flatint ' -o ' flatint ' -r 1 -a 1\n'];
        fprintf(fid_run,command);
        %command=['imageMath.py -e=''abs(b)*exp(1.0*J*arg(b*a*conj(b)*c))'' -o ' flatslc ' -t cfloat --a=' primary ' --b=' secondary ' --c=' flatint '\n'];
        command=['imageMath.py -e=''b*a*conj(b)*c/abs(c)**2'' -o ' flatslc ' -t cfloat --a=' primary ' --b=' secondary ' --c=' flatint '\n'];

        fprintf(fid_run,command);
%         if(i>1)
%             fprintf(fid_run,'rm -r %s\n',intdir);
%         end
    end
end 
fclose(fid_run);
