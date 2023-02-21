function local_stiffness = new_localstiff(node, Po ,Lege_ind, a)



dim_elem =size(Lege_ind,1); % number of basis for each element.

local_stiffness = zeros(dim_elem,dim_elem);

[w_x,x] = quad_GL(ceil((Po+1)*0.5));   [w_y,y] = quad_GJ1(ceil((Po+1)*0.5));

quad_x = kron(x,ones(size(w_y,1),1)); quad_y = kron(ones(size(w_x,1),1),y);

 
 shiftpoints = [(1+quad_x).*(1-quad_y).*0.5-1, quad_y ];

 ref_points = 0.5.*shiftpoints+0.5;  

% For each triangle

    %B is the affine matrix
    
    B = 0.5*[node(2,:)-node(1,:); node(3,:)-node(1,:)];
    
    De_tri = abs(det(B));   inv_B = inv(B);
    
    weights = kron(w_x,w_y) ; weights = weights.*0.5.*De_tri;
              
    phy = reference_to_physical_t3 (node([1,2,3],:)', size(ref_points,1), ref_points' );
    
    P_Qpoints  = phy';       
    
    a_val = a(P_Qpoints); 

    % data for quadrature
  

stiff = zeros(dim_elem,dim_elem);

for i = 1:dim_elem
   
    %for j = 1:dim_elem
    
    %%symetric term
    
    for j=i:dim_elem
    
        % first term {a grad(u) \cdot grad(v)} is symetric

        
        %% gradient over reference element
        
        grad1 = GradSimplex2DP(quad_x,quad_y,Lege_ind(i,1),Lege_ind(i,2),inv_B);
        
        grad2 = GradSimplex2DP(quad_x,quad_y,Lege_ind(j,1),Lege_ind(j,2),inv_B);
        
        
        
        grad = [grad1(:,1).*grad2(:,1) , grad1(:,1).* grad2(:,2),grad1(:,2).*grad2(:,1) , grad1(:,2).* grad2(:,2)];
        
        t = sum((grad).*a_val,2);
        
       
        stiff(j,i) =dot(t,weights);
        
    end
end

%%symetric term



local_stiffness = local_stiffness + stiff + (tril(stiff,-1))';

end

