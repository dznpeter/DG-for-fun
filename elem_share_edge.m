 function edge2elem = elem_share_edge(input_edge,totalEdge,elemtotaledge)
 
  T = totalEdge - kron(input_edge,ones(size(totalEdge,1),1))  ;
 
  T = sum(abs(T),2); 
  
 [index  , ~ , ~] = find(T==0);
 
 
 edge2elem = NaN(size(index,1),1);
 
 for i = 1:size(index,1)
 
      t = elemtotaledge -index(i);
      
      [k  , ~ , ~] = find(min(t,0) ==0);
      
      edge2elem(i) = k(1);
     
 end
 
 end