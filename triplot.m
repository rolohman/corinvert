function triplot(data,dn)
nd=length(dn);
cids=[];
for i=1:nd
    cids=[cids (nd)*(i-1)+[i:nd-1]];
end

jnk=nan(nd);
jnk(cids)=data;
pcolor(dn,dn,jnk'),shading flat,set(gca,'ydir','reverse');
