function geocoder(file,pol,resdir,bbox,dem,geodir)
disp(['geocoding '  file ' ' pol])
decide_ints_stack
chdir(resdir)

command=['imageMath.py -e ''a'' -o tmp --a=''' file ';' num2str(newnx) ';float;1;BSQ'''];
system(command);
copyfile(file,'tmp');

command=['geocodeIsce.py -f tmp -b ''' num2str(bbox) ''' -d ' dem ' -m ../master -s ../master -r ' num2str(rlooks) ' -a ' num2str(alooks)];
system(command);

movefile('tmp.geo',[file '.geo']);
movefile('tmp.geo.vrt',[file '.geo.vrt']);
movefile('tmp.geo.xml',[file '.geo.xml']);
command=['fixImageXml.py -i ' file '.geo -b'];
system(command);

chdir('..');