function val_p  = shift_leg_derivative(x,m,h,order,k)


tol =10^(-15);

y = (x-m)./h;  %shift to the interval [-1,1]%


test = abs(y)-1;

[ia,~,~] = find(test >= 0);

if isempty(ia) ==0
    
    
   y(ia,:)=(1-tol).*(y(ia,:)./abs(y(ia,:)));
    
end


% Tolerance for Legendre function


if order >= 0 && order <= k-1 
   
    val_p = 0.*y;
    
else
    
P = JacobiP(y,k,k,order-k).*sqrt(gamma(order+k+1)./gamma(order-k+1)); 

val_p = h.^(-0.5-k).*P;
    


end


end