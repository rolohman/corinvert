function copyxml(in,out)
copyfile([in '.xml'],[out '.xml']);
copyfile([in '.vrt'],[out '.vrt']);
command=['fixImageXml.py -i ' out ' -f'];
system(command);