function z = diffusive(x)


if size(x,2)~= 2
    
   error('input should be 2 dimensional points')
    
end

%{

eps=0.1;

z = [eps+x(:,1),   x(:,1).*x(:,2)   ,  x(:,1).*x(:,2)   ,   eps+x(:,2) ];

%}

z = [ ones(size(x,1),1),   zeros(size(x,1),1)   , zeros(size(x,1),1)  ,   ones(size(x,1),1) ];

end