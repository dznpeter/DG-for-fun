%plot the functions

clear all;

x=-1:0.05:1;

[xx, yy]=meshgrid(x,x);

zz=u_true([xx(:),yy(:)]);

surf(xx,yy,reshape(zz, size(xx,1),[]))
