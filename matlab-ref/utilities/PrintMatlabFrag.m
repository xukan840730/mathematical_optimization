function [] = PrintMatlabFrag(ncols,vncols,fontsize,fname)

%set horizontal figure extent to fit in 1/2 or whole typical page width
hsize=8.7.*ncols;

vsize=8.7.*vncols; %set vertical figure extent, in units of 'column width'

set(gcf,'units','centimeters');
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2),hsize,vsize]);

%Set all fontsizes
set(findall(gcf,'-property','FontSize'),'FontSize',fontsize);

%print to .eps + .tex files
matlabfrag(fname);

end