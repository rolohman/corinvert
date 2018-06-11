function cr=est_cr(cres,synth,nd,cidl)

crs              = zeros(nd);
crs(cidl)        = cres;
crs              = crs+crs';
crs(crs==0)      = NaN;

   
if(isempty(synth))
    cerr=0.3083*(1-crs)-.2*(1-crs).^2;
    wgts=1./cerr;
else
    wgts             = zeros(nd);
    wgts(cidl)       = synth;  %0 to 1
    wgts             = wgts+wgts';
    wgts(wgts==0)    = NaN;
    wgts             = exp((wgts-1)*5); %scalar makes near 0.1 at synth=0.5;
end
cr=mymed(crs,wgts);
