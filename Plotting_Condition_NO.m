
%% plotting errors for the DG for p refinement

Polynomial_degree = 1:10;    penalty=10;

NO_elem = [32];

Condition_NO_vect =NaN(length(Polynomial_degree),1); Polynomial_vect = NaN(length(Polynomial_degree),1);

Mesh_type = 'triangle'; 

figure;

for i=1 :length(Polynomial_degree)
    
 
load(['Error ' num2str(NO_elem) ' ' Mesh_type ' Elements penalty ' num2str(penalty) ' P' num2str(Polynomial_degree(i)) ' basis.mat'])

        
Condition_NO_vect(i)=Condition_NO;      
       

Polynomial_vect(i) = i;

end



log_Condition_NO_vect = log(abs(Condition_NO_vect)); 

slope_P1 = abs((log_Condition_NO_vect(2:end)-log_Condition_NO_vect(1:end-1))./(log(Polynomial_vect(2:end))-log(Polynomial_vect(1:end-1))));


%slope_Poly(k) = max(mean(slope_P1(end-1:end)), max(slope_P1));

slope_Poly = mean(slope_P1(end-3:end));

%slope_Poly = max(slope_P1(end));






loglog(Polynomial_vect,Condition_NO_vect,'r-s','LineWidth',2,'MarkerSize',10);

legend(['DG condition No order  ' num2str(slope_Poly)]...
      ,'Location','NorthWest')
  


xlabel('p','FontSize',18);

        
ylabel('Condition Number','FontSize',20);      
       



set(gca,'FontSize',20)

%title('P refinement','FontSize',20)




