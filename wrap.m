function  [index, swap] = wrap(elem_oneface,tri1)


edge = elem_oneface';

%index = 0;  swap = 0;

edge_tri = [tri1, tri1([2,3,1])];


test1 = abs((edge_tri(:,1)-edge(1))) +abs((edge_tri(:,2)-edge(2)));

index = find(test1==0);

if isempty(index) == 0        
    
    swap = num2str(0);
else
    
    swap = num2str(1);
    
    test1 = abs((edge_tri(:,1)-edge(2))) +abs((edge_tri(:,2)-edge(1)));
    
    index = find(test1==0);
    
end
   
  
index = num2str(index);









