%% test for Jacobi polynomial 

%{

n=5;

[weights, ref_points] = quad_GL(n^2);

alpha = 8;  beta=0;


P1 = Deriv_JacobiP(-1,alpha,beta,n,0);

P2 = Deriv_JacobiP(1,alpha,beta,n,0);

P2 - P1

 v = Deriv_JacobiP(ref_points,alpha,beta,n,1);
 
 dot(v,weights)
 
 %}

test_B = eye(2);

a = -0.2;

b = -0.4;

i = 19; 

j = 4;


grad = GradSimplex2DP(a,b,i,j, test_B);


[dmodedr, dmodeds] = test_GradSimplex2DP(a,b,i,j);

grad - [dmodedr, dmodeds]




