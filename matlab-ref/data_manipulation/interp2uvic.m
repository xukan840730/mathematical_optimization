function[arrayout] = interp2uvic(arrayin)

%routine interpolates an arbitrary regular 2D array to the UVic grid.

dim = size(arrayin);

xi=linspace(1,dim(2),100);
yi=linspace(1,dim(1),100);
xo=1:1:dim(2);
yo=1:1:dim(1);

[X,Y] = meshgrid(xo,yo);
[XI,YI] = meshgrid(xi,yi);  
arrayout=interp2(X,Y,arrayin,XI,YI,'*spline');