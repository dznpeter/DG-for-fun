function [quad_x, quad_y] = quadpreference(ref_Qpoints, index, swap)

% reference triangel is (-1,-1)--(1,-1)--(-1,1), 
% Edge 1 is (-1,-1)--(1,-1)
% Edge 2 is (1,-1)--(-1,1)
% Edge 2 is (-1,1)--(-1,-1)



if    index == '1'
        
    quad_x = ref_Qpoints;     
    
    quad_y = -ones(size(ref_Qpoints)); 
    
    
end


if  index == '2'
        
     quad_x = -ref_Qpoints; 
     
     quad_y = ref_Qpoints; 
     
end

     
if   index == '3'       
        
      
     quad_x = -ones(size(ref_Qpoints));     
    
     quad_y = -ref_Qpoints;    
     
    

end
        

if swap == '1'
    
    quad_x = quad_x(end:-1:1);
    
    quad_y = quad_y(end:-1:1);
end


end