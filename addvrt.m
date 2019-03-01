function addvrt(file,nx)
command=['imageMath.py -e ''a'' -o tmp --a=''' file ';' num2str(nx) ';float;1;BSQ'''];
system(command);
movefile('tmp.vrt',[file '.vrt']);
movefile('tmp.xml',[file '.xml']);
command=['fixImageXml.py -i ' file ' -f'];
system(command);
