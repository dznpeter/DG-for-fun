function latextable_compare_error(filename)
%% LATEXTABLE output the matrix in table formate in latex.
%
% latextable(A,'A.tex');
% latextable(A);
% 





if (nargin <1), filename='table.tex'; end %default file name
fid = fopen(filename,'wt');

%% Latex
fprintf(fid,'\\documentclass[11pt]{amsart} \n'); 
fprintf(fid,'\\usepackage{geometry} \n'); 
fprintf(fid,'\\geometry{letterpaper} \n'); 
fprintf(fid,'\\usepackage{graphicx}  \n'); 
fprintf(fid,'\\usepackage{amssymb}  \n'); 
fprintf(fid,'\\usepackage{epstopdf}  \n'); 
fprintf(fid,'\\DeclareGraphicsRule{.tif}{png}{.png}{`convert #1 `dirname #1`/`basename #1 .tif`.png} \n\n \n'); 
fprintf(fid,' \\title{DG VS VEM} \n\n \n'); 
fprintf(fid,' \\author{The Author}\n\n\n '); 
fprintf(fid,' \\begin{document}\n'); 
fprintf(fid,' \\maketitle\n\n\n '); 



%% k is the order of polynomial basis

for k = 1 : 4


%% load the data of error

penalty=[1,1.5,2,3,5];    

NO_elem = [12,48,192,768,3072,12288,49152 ];


n = 2+length(penalty);    
    
data = cell(1,n);  

data{1} = '\\text{Cell} '; data{2} = '\\text{Dofs}';

data_H1 = data;

for i = 3:n

data{i} =['\\text{error $C_\\sigma=$}' num2str(penalty(i-2)) ]; 

data_H1{i} =['\\text{error $C_\\sigma=$}' num2str(penalty(i-2)) ]; 

end




m = length(NO_elem);


A = NaN(m,n);

B = NaN(m,n);

for i=1 :length(NO_elem)
    
    

load(['Error ' num2str(NO_elem(i)) ' rectangle Elements penalty ' num2str(penalty(1)) ' P' num2str(k) ' basis.mat'])



A(i, 1) =NO_elem(i) ;  A(i, 2) =dim_FEM ;   A(i,3) = L2_err;


B(i, 1) =NO_elem(i) ;  B(i, 2) =dim_FEM ;   B(i,3) = H1_err;



for j=4:n
    
    load(['Error ' num2str(NO_elem(i)) ' rectangle Elements penalty ' num2str(penalty(j-2)) ' P' num2str(k) ' basis.mat'])

    
    A(i,j) = L2_err;
    
    B(i,j) = H1_err;
    
end


end






%% Latex file for Table for L2 error

fprintf(fid,'\\begin{table}[!htb] \n');

fprintf(fid,'\\begin{center}\n');

fprintf(fid,'\\begin{tabular}{');
for j=1:n
    fprintf(fid,'|c');
end
fprintf(fid,'|} \n \\hline  \n');


for j=1:n
        fprintf(fid,data{j});
        if j<n
            fprintf(fid,' & ');
        end
end
    fprintf(fid,'\\\\ \n \\hline   \\hline \n');
    
    
   

for i=1:m
    for j=1:n
        
        if j <=2
        
        fprintf(fid,'%g',A(i,j));
        
        else
            
        fprintf(fid,'%e',A(i,j));    
        
        end
        if j<n
            fprintf(fid,' & ');
        end
    end
    fprintf(fid,'\\\\ \n \\hline \n');
end

fprintf(fid,'\\end{tabular}  \n ');

fprintf(fid,'\\end{center}\n');

fprintf(fid,['\\caption{DG method Polynomial Order $k=$' num2str(k) ', $\\bold{L_2}$ norm error with different penalty}\n ']);

fprintf(fid,'\\end{table} \n  \n');





%% Latex file for Table for H1 error

fprintf(fid,'\\begin{table}[!htb] \n');

fprintf(fid,'\\begin{center}\n');

fprintf(fid,'\\begin{tabular}{');
for j=1:n
    fprintf(fid,'|c');
end
fprintf(fid,'|} \n \\hline  \n');


for j=1:n
        fprintf(fid,data_H1{j});
        if j<n
            fprintf(fid,' & ');
        end
end
    fprintf(fid,'\\\\ \n \\hline   \\hline \n');
    
    
   

for i=1:m
    for j=1:n
        
        if j <=2
        
        fprintf(fid,'%g',B(i,j));
        
        else
            
        fprintf(fid,'%e',B(i,j));    
        
        end
        if j<n
            fprintf(fid,' & ');
        end
    end
    fprintf(fid,'\\\\ \n \\hline \n');
end

fprintf(fid,'\\end{tabular}  \n ');

fprintf(fid,'\\end{center}\n');

fprintf(fid,['\\caption{DG method Polynomial Order $k=$' num2str(k) ', $\\bold{H_1}$ semi-norm error with different penalty}\n ']);

fprintf(fid,'\\end{table} \n  \n');


end

fprintf(fid,'\\end{document} \n  \n');
fclose(fid);

%% TODO: add latexerrtable