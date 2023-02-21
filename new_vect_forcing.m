function z = new_vect_forcing(node,  Po, Lege_ind, f)

dim_elem =size(Lege_ind,1); % number of basis for each element.

z = zeros(dim_elem,1);  % initialize the local vector


%[w_x,x] = Golub_Welsch(ceil(Po*0.5));

[w_x,x] =quad_GL(ceil((Po+3)*0.5));

%[w_y,y] = JacobiGQ(1,0,ceil(Po*0.5));

[w_y,y] = quad_GJ1(ceil((Po+3)*0.5));
 
%[w_y,y] = Golub_Welsch(ceil((Po+1)*0.5));

quad_x = kron(x,ones(size(w_y,1),1)); quad_y = kron(ones(size(w_x,1),1),y);
 
 
 %quadpoins = [quad_x,quad_y];
 
weights = kron(w_x,w_y) ;

shiftpoints = [(1+quad_x).*(1-quad_y).*0.5-1, quad_y ];
 
 ref_points = 0.5.*shiftpoints+0.5;  


 % For each triangle


    
    %B is the affine matrix
    
    B = 0.5*[node(2,:)-node(1,:); node(3,:)-node(1,:)];
              
        De_tri = abs(det(B)); weights = weights.*0.5.*De_tri;

    
    phy = reference_to_physical_t3 (node([1,2,3],:)', size(ref_points,1), ref_points' );
    
    P_Qpoints  = phy';
    
 
     
    % data for quadrature
    
    %%% numerical intergration by Duffy Transfomation
     
       
    f_val = f(P_Qpoints); 

first = zeros(dim_elem,1);
    
for j = 1:dim_elem
   
   
    
        % first term {fv} 
         
        %t = f_val.*tensor_leg(P_Qpoints,m,h,Lege_ind(j,:)).*Duffy_y;
        
        
        t = f_val.*Simplex2DP(quad_x,quad_y,Lege_ind(j,1),Lege_ind(j,2));
        
        
        first(j) = dot((t),weights);
      
    
end

z = z + first;









end