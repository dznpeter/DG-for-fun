clear all; close all;

%% different elements we use different order of basis

NO_Elements = [32,128,512,2048,8192];


color = {'b','r','g','k'};

L2=figure;

slope_Rect = NaN(1,4);

order = 4;

for k = 1:order


%% plotting errors for the L2 norm for p refinement

%%P basis


Polynomial_degree = k;    penalty=10;

NO_elem = NO_Elements;   

err_P = NaN(length(NO_elem),1); dof_P = NaN(length(NO_elem),1);


for i=1 :length(NO_elem)

load(['Error ' num2str(NO_elem(i)) ' triangle Elements penalty ' num2str(penalty) ' P' num2str(Polynomial_degree) ' basis.mat'])
    
    
err_P(i,1)=L2_err;         dof_P(i) = dim_FEM;

       

end



%slope for convergence line

logerr_P1 = abs(log(err_P(:,1))); 

slope_P1 = abs((logerr_P1(2:end)-logerr_P1(1:end-1))./(log(dof_P(2:end).^(1./2))-log(dof_P(1:end-1).^(1./2))));



slope_Rect(k) =  mean(slope_P1(4:end));



loglog(dof_P,err_P(:,1),['-s' num2str(color{k})],'LineWidth',1.5,'MarkerSize',10);

hold on;

end




legend(['P1 slope ' num2str(slope_Rect(1))]...
      ,['P2 slope ' num2str(slope_Rect(2))]...
      ,['P3 slope ' num2str(slope_Rect(3))]...
      ,['P4 slope ' num2str(slope_Rect(4))]...
      ,'Location','SouthWest')


xlabel('$\sqrt{\rm DoFs}$','FontSize',20,'Interpreter','latex');
        
  
ylabel('$||u-u_{h}||_{L^2(\Omega)}$','FontSize',20,'Interpreter','latex');
 




set(gca,'FontSize',20)



 %% plotting errors for the DG norm for P refinement
 
 %%P basis
 
 H1=figure;


for k = 1:order


%% plotting errors for the H1 norm for h refinement

%%P basis


Polynomial_degree = k;    penalty=10;

NO_elem = NO_Elements;   

err_P = NaN(length(NO_elem),1); dof_P = NaN(length(NO_elem),1);


for i=1 :length(NO_elem)
       

    
load(['Error ' num2str(NO_elem(i)) ' triangle Elements penalty ' num2str(penalty) ' P' num2str(Polynomial_degree) ' basis.mat'])

err_P(i,1)=DG_err;         dof_P(i) = dim_FEM;
        


end

%slope for convergence

logerr_P1 = abs(log(err_P(:,1))); 

slope_P1 = abs((logerr_P1(2:end)-logerr_P1(1:end-1))./(log(dof_P(2:end).^(1./2))-log(dof_P(1:end-1).^(1./2))));

slope_Rect(k) = mean(slope_P1(4:end));



loglog(dof_P,err_P(:,1),['-s' num2str(color{k})],'LineWidth',1.5,'MarkerSize',10);

hold on;


end




legend(['P1 slope ' num2str(slope_Rect(1))]...
      ,['P2 slope ' num2str(slope_Rect(2))]...
      ,['P3 slope ' num2str(slope_Rect(3))]...
      ,['P4 slope ' num2str(slope_Rect(4))]...
      ,'Location','SouthWest')




xlabel('$\sqrt{\rm DoFs}$','FontSize',20,'Interpreter','latex');
        
  
ylabel('$|||u-u_{h}|||_{\rm DG}$','FontSize',20,'Interpreter','latex');
 



set(gca,'FontSize',20)

