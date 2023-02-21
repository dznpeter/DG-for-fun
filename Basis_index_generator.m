function  Lege_ind = Basis_index_generator(Polydegree,Dimension)

if Dimension ==2
   
               

% P basis

t=0:Polydegree;  t=t';

temp = ones(Polydegree+1,1);

Lege_ind = [kron(t,temp) , kron(temp,t)];

% delete the index 

index = find( sum(Lege_ind,2) > Polydegree );

Lege_ind(index,:) = [];

end





if Dimension ==3
   

% P basis

t=0:Polydegree;  t=t';

temp = ones(Polydegree+1,1);

Lege_ind = [  kron(kron(t,temp),temp) , kron(kron(temp,t),temp), kron(kron(temp,temp),t)];


% delete the index 

index = find( sum(Lege_ind,2) > Polydegree );

Lege_ind(index,:) = [];


end


end