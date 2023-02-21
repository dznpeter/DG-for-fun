function grad = GradSimplex2DP(a,b,i,j, inv_B)


%% inverse B is used to map from physical tri to reference triangle


%val_p  = Deriv_JacobiP(x,order,k,alpha,beta)

h1 = JacobiP(a,0,0,i);  dh1 = Deriv_JacobiP(a,0,0,i,1);

h2 = JacobiP(b,2*i+1,0,j);  dh2 = Deriv_JacobiP(b,2*i+1,0,j,1);

%P = sqrt(2.0)*h1.*h2.*(1-b).^i;



grad1 = sqrt(2.0)*dh1.*h2.*(1-b).^i.*(2./(1-b)) ;

grad2 = sqrt(2.0)*h1.*dh2.*(1-b).^i -sqrt(2.0)*h1.*h2.*(1-b).^(i-1).*i...
       +sqrt(2.0)*dh1.*h2.*(1-b).^i.*((1+a)./(1-b)) ;

%%

        inv_gradu1 = inv_B(1,1).*grad1 + inv_B(1,2).*grad2;
        
        inv_gradu2 = inv_B(2,1).*grad1 + inv_B(2,2).*grad2;
        

        grad = [inv_gradu1, inv_gradu2];
        
     