% the Golub-Welsch algorithm with Legendre Polynomial


function [w,x] = Golub_Welsch(n)

% to transfer the quadrature for interval [a,b] to [-1,1]

%off-diagonal entries of Legendre's L

b=1./(2*sqrt(1-(2*(1:n)).^(-2)));

%the symmetric tridiagonal matrix L of size n+1

L=diag(b,1)+diag(b,-1);

%V is eignevector and X is eigenvalue
[V,X]=eig(L);
%Gauss nodes xj are eigenvalues of L
x=diag(X);
[x,j]=sort(x);
%the weights w are positive
w=2*V(1,j).^2;
w=w';




end
