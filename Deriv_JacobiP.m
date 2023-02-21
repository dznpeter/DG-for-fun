function val_p  = Deriv_JacobiP(x,alpha,beta,order,k)



if order >= 0 && order <= k-1 
   
    val_p = 0.*x;
    
else
    
    
c1 = gamma(order+1).*gamma(order+k+alpha+beta+1);  

c2 = gamma(order-k+1).*gamma(order+alpha+beta+1); 
    
val_p = JacobiP(x,alpha+k,beta+k,order-k).*sqrt(c1/c2); 

end


end