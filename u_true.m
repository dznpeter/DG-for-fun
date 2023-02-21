function z = u_true(x)


if size(x,2)~= 2
    
   error('input should be 2 dimensional points')
    
end


 %Case 1

  z = sin(pi.*x(:,1)).*sin(pi.*x(:,2));


% %Case 2
% 
% 
% r2 = (x(:,1)-0.5).^2+(x(:,2)-0.5).^2;
% 
% z = exp(-r2);

% %Case 3

%z = sin(pi.*x(:,1)).*sin(2*pi.*x(:,2));

end