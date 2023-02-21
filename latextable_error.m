function latextable_error(filename)
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


%% Latex for figures

fprintf(fid,'\\begin{figure} [!htb] \n');

fprintf(fid,'\\begin{center}\n');

fprintf(fid,'\\begin{tabular}{');

for j=1:2
    fprintf(fid,'c');
end
fprintf(fid,'}   \n');

fprintf(fid,'\\includegraphics[width=0.48 \\linewidth]{Poisson_L_shpae_L2_norm_error_h_refine} &  \\includegraphics[width=0.48 \\linewidth]{Poisson_L_shpae_H1_semi_norm_error_h_refine}\n');


fprintf(fid,'\\end{tabular}  \n ');

fprintf(fid,'\\end{center}\n');

fprintf(fid,'\\caption{$L_2$ norm (left) and $H_1$ semi-norm (right) error under $h$  refinement for fixed $p$} \n ');

fprintf(fid,'\\end{figure} \n  \n');


%% k is the order of polynomial basis

for k = 1 : 4

n = 6;    
    
data = cell(1,n);

data{1} = '\\text{NO of Elements} '; data{2} = '\\text{Dofs}';

data{3} ='$L_2$\\text{ norm error}'; data{4} ='$L_2$\\text{ rate } $(h)$';    

data{5} = '$H_1$\\text{ semi-norm error}'; data{6} ='$H_1$\\text{ rate} $(h)$';  


%% load the data of error

penalty=10;  NO_elem = [12,48,192,768,3072,12288,49152,196608 ];


m = length(NO_elem);


A = NaN(m,n);



for i=1 :length(NO_elem)

load(['Error ' num2str(NO_elem(i)) ' rectangle Elements penalty ' num2str(penalty) ' P' num2str(k) ' basis.mat'])



A(i, 1) =NO_elem(i) ;  A(i, 2) =dim_FEM ;


if i  ==  1
   
A(i, 3) = L2_err;      

A(i, 5) = H1_err;     
    
    
    
else


  
A(i, 3) = L2_err;      A(i,4) = 2.*(abs(log(A(i, 3)))-abs(log(A(i-1, 3))))./(log(NO_elem(i))-log(NO_elem(i-1)));

A(i, 5) = H1_err;      A(i,6) = 2.*(abs(log(A(i, 5)))-abs(log(A(i-1, 5))))./(log(NO_elem(i))-log(NO_elem(i-1)));



end
end






%% Latex file for Table

fprintf(fid,'\\begin{table}[!htb]  \n');

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
        
        if j == 3 || j==5 
        
        fprintf(fid,'%e',A(i,j));
        
        else
            
        fprintf(fid,'%g',A(i,j));    
        
        end
        if j<n
            fprintf(fid,' & ');
        end
    end
    fprintf(fid,'\\\\ \n \\hline \n');
end

fprintf(fid,'\\end{tabular}  \n ');

fprintf(fid,'\\end{center}\n');

fprintf(fid,['\\caption{DG method Polynomial Order $k=$' num2str(k) '}\n ']);

fprintf(fid,'\\end{table} \n  \n');
end

fprintf(fid,'\\end{document} \n  \n');
fclose(fid);

%% TODO: add latexerrtable