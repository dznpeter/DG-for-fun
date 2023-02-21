function z= forcing(x)


if size(x,2)~= 2
    
   error('input should be 2 dimensional points')
    
end

% Diffussion part



% eps=0.1;
% 
% t1 = -pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2)) + pi.^2.*(eps+x(:,1)).*sin(pi.*x(:,1)).*sin(pi.*x(:,2));
% 
% t2 = -pi.*x(:,1).*cos(pi.*x(:,1)).*sin(pi.*x(:,2)) - pi.^2.*x(:,1).*x(:,2).*cos(pi.*x(:,1)).*cos(pi.*x(:,2));
% 
% t3 = -pi.*x(:,2).*cos(pi.*x(:,2)).*sin(pi.*x(:,1)) - pi.^2.*x(:,1).*x(:,2).*cos(pi.*x(:,1)).*cos(pi.*x(:,2));
% 
% t4 = -pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1)) + pi.^2.*(eps+x(:,2)).*sin(pi.*x(:,1)).*sin(pi.*x(:,2));
% 
% z = (t1+t2+t3+t4);

 % case1
 
 
  z = 2.*pi.^2.*sin(pi.*x(:,1)).*sin(pi.*x(:,2));



% % case 2
% 
% r2 = (x(:,1)-0.5).^2+(x(:,2)-0.5).^2;
% 
% z =4.*(1-r2).*exp(-r2) ;
 
 
 % %Case 3

%z = [pi*cos(pi.*x(:,1)).*sin(0.5*pi.*x(:,2)), 0.5*pi*sin(pi.*x(:,1)).*cos(0.5*pi.*x(:,2))];

%z = 5.*pi.^2.*sin(pi.*x(:,1)).*sin(2*pi.*x(:,2));

end