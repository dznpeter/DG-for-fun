function z = grad_u_true(x)


if size(x,2)~= 2
    
   error('input should be 2 dimensional points')
    
end

 %Case 1

 z = [pi*cos(pi.*x(:,1)).*sin(pi.*x(:,2)), pi*sin(pi.*x(:,1)).*cos(pi.*x(:,2))];


% % %Case 2
% 
 %r2 = (x(:,1)-0.5).^2+(x(:,2)-0.5).^2;
% 
 %z = [exp(-r2).*(-2.*(x(:,1)-0.5)) , exp(-r2).*(-2.*(x(:,2)-0.5))];

% %Case 3

%z = [pi*cos(pi.*x(:,1)).*sin(2*pi.*x(:,2)), 2*pi*sin(pi.*x(:,1)).*cos(2*pi.*x(:,2))];


end