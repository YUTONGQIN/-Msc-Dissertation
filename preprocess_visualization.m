%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%import normalized gene counts data
%the ensembl ID should be importred as a string vector
gene_counts = MayoRNAseqRNAseqTCXgeneCountsnormalized;
gene_counts = table2cell(gene_counts);

%transpose data matrix
gene_counts = cell2table(gene_counts');

%check the table contents
gene_counts(1:5,1:5)

%create a string vecoter containing the row names of this table
rownames = cellstr(string(gene_counts{:,1}));
rownames(1) = [];

%create a string vecoter containing the column names of this table 
varnames = cellstr(gene_counts{1,:});

%delete the first row of this table
gene_counts(1,:) = [];

%transfer each column of table from string to double
genecounts_matrix =[];
for i=1:size(gene_counts,2)
  col = str2double(gene_counts{:,i});
  genecounts_matrix = [genecounts_matrix col];
  disp(i)
end

%remove columns with all zero elements
genecounts_matrix_old = genecounts_matrix;
genecounts_matrix(:,~any(genecounts_matrix,1)) = [];

%transfer the matrix b back to table and add both row and column names 
gene_counts = array2table(genecounts_matrix);
newvarnames = varnames(any(genecounts_matrix_old,1));

gene_counts.Properties.RowNames = rownames;
gene_counts.Properties.VariableNames = newvarnames;

%save 'gene_counts' as a txt file so that we can use it whenever we need
%writetable(gene_counts,'C:\Work Files\project\data\gene_counts.txt')






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%import data containing information about covariates  
covariates = MayoRNAseqRNAseqTCXcovariates;

%find the row name of 'covariates'
rownamesB = cellstr(string(covariates{:,1}));
covariates.Properties.RowNames = rownamesB;

%check the table
head(covariates)

%change the column name ID to enxemble_id so that we can merge tables based on a common column
covariates.Properties.VariableNames{'ID'} = 'ensembl_id'; 

%save 'covariates' as a txt file
%writetable(covariates,'C:\Work Files\project\data\covariates.txt')






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%only keep variables required
covariates = covariates(:,[1 4 5 6 7]);

%merge 'gene_counts' and 'covariates'
merged_table = join(gene_counts,covariates,'Keys','ensembl_id');
merged_table(:,1)=[];

%delete rows with missing values (these missing values come from 'covariates')
merged_table = merged_table(~any(ismissing(merged_table),2),:);

%only keep AD and Control samples
index_AD = find(merged_table{:,48230} =='AD');
index_Ctrl = find(merged_table{:,48230} =='Control');
index = [index_AD;index_Ctrl];
newmerged_table = merged_table(index,:);

%check table
newmerged_table(1:100,48227:48232)

%record the rows and columns name of the new merged table
a = newmerged_table.Properties.VariableNames;
b = newmerged_table.Properties.RowNames;

%regress three covariates(sex, RNA quality and age) out
tb_aftctl = [];
for i = 1:48228
     table = newmerged_table(:,[48229 48231 48232 i]);
     mdl = fitlm(table);
     tb_aftctl = [tb_aftctl mdl.Residuals.Raw];
     disp(i)
end

%add sample lables (ADs are labelled by 1 and Controls 0)
Diagnosis = double(newmerged_table{:,48230}=="AD");
tb_aftctl = [tb_aftctl Diagnosis];

%construct original matrix
raw_data = table2array(newmerged_table(:,1:48228));
tb_aftctl2 = [raw_data Diagnosis];

%check the proportion of AD samples
sum(Diagnosis==1)
[Diagnosis reshape(1:length(Diagnosis),160,1)]

%delete columns with all zero elements
index_delete = ~any(tb_aftctl,1);
tb_aftctl(:,index_delete) = [];

%save 'tb_aftctl' as a csv file
%dlmwrite('C:\Work Files\project\data\tb_aftctl.csv',tb_aftctl);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%draw MDS plots using different distances
%calculate dissimilarity matrix
D = pdist(tb_aftctl(:,1:47201),'cityblock');
D = pdist(tb_aftctl(:,1:47201),'euclidean');
D = pdist(tb_aftctl(:,1:47201),'correlation');

%non-metric MDS (force dimensionality to 2D)
Y = cmdscale(D, 2)';

%display results
figure;
gca = plot(Y(1, 1:82), Y(2, 1:82), 'b.','markersize',20);
hold on;
grid on;
plot(Y(1, 83:160), Y(2, 83:160), 'r.','markersize',20);
legend({'AD' 'Control'},...
       'Location','NorthEast');
hold off;
 
%non-metric MDS (force dimensionality to 3D)
Y = cmdscale(D, 3)';

%display results
figure;
gca2 = plot3(Y(1, 1:82), Y(2, 1:82),Y(3, 1:82), 'b.','markersize',20);
hold on;
plot3(Y(1, 83:160), Y(2, 83:160),Y(3, 83:160), 'r.','markersize',20);
grid on
legend({'AD' 'Control'},...
       'Location','NorthEast');
hold off;   