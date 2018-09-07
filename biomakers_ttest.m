%replace the original class lablels (1 and 0) by new ones (1 and -1)
label_set = tb_aftctl(:,47202);
label_set = sign(label_set-0.5);

%generate five-fold cross-Validation indices of the outer loop
Indices2 = crossvalind('Kfold',label_set, 5);

output = [];
LOOE = cell(1,5);
LOOE2 = cell(1,5);

Avg_log_LS = cell(1,5);
Features_Hits = cell(1,5);
Rank_in_ttest = cell(1,5);

for ii = 1:5
 
%four folds are used as training and validation data
train_valid_set = tb_aftctl(Indices2~=ii,1:47201);
train_valid_label = label_set(Indices2~=ii);

%one fold left out is regarded as test dataset
test_set = tb_aftctl(Indices2==ii,1:47201);
test_label = label_set(Indices2==ii);

%generate five-fold cross-Validation indices of the inner loop
Indices = crossvalind('Kfold',train_valid_label, 5);

recorded_labels = cell(1,5);
recorded_LS = cell(1,5);

%carry out ttest in pre-processing to remove the most uninformative genes
[h,pvalue,ci,stat] = ttest2(train_valid_set(train_valid_label==1,:),train_valid_set(train_valid_label==-1,:),'Vartype','unequal'); 
[sorted_pvalue,index1] = sort(pvalue);
 
%select the first 100 features with smaller p-values
IDX_t = index1(1:100); 

for i = 1:5

%four folds are used as training data
train_set = train_valid_set(Indices~=i,IDX_t);
train_label = train_valid_label(Indices~=i);

%one fold left out is regarded as validation dataset
valid_set = train_valid_set(Indices==i,IDX_t);
valid_label = train_valid_label(Indices==i);

%find the number of dimensions (D) of the training data
[n,D] = size(train_set);

%set the log initial values of D+1 hyperparameters (D length scales and a signal
%variance)
%set each length scale to be the average of pairwise Euclidean distances 
ell = repmat(mean(mean(pdist2(train_set,train_set))),1,D);      

%set the signal variance to be 1
sf = 1;
loghyper = log([ell sf]);
 
%maximize the log marginal likelihood w.r.t hyperparameters by carrying out
%conjugate gradient optimization (-50 tells minimize to evaluate the function at most 50 times)
[loghyper binaryLaplaceGP(loghyper, 'covSEard', 'cumGauss', train_set, train_label)];
[newloghyper logmarglik] = minimize(loghyper, 'binaryLaplaceGP', -50, 'covSEard', 'cumGauss', train_set, train_label);
[newloghyper logmarglik(end)];
     
%sort the log length scales in ascending order (equivalent to sort the ARD values in descending order)
[sorted_nloghyper,index2] = sort(newloghyper(1:D));

%rank features based on the optimal log length scales (the optimal ARD values)
data_set1 = train_set(:,index2);
data_set_label1 = train_label;

%perform forward selection to calculate Leave-One-Out errors
E = forwardFS_LOO(data_set1,data_set_label1);

%find top-ranked features achieving the minimal LOO error 
min_loc = find(E==min(E))-1;  
origin = 1:47201;
n_origin = origin(index1);
nn_origin = n_origin(1:100);
nnn_origin = nn_origin(index2);
selected_features =  nnn_origin(1:min_loc);  
 
%record the selected features and their log length scales
recorded_labels{i} = selected_features; 
recorded_LS{i} = sorted_nloghyper(1:min_loc);

end
  





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%combine features selected in each fold in the inner loop
feature_combination = [recorded_labels{1} recorded_labels{2} recorded_labels{3} recorded_labels{4} recorded_labels{5}];

%combine the log length scales of features selected in each fold in the inner loop
LS_combination = [recorded_LS{1} recorded_LS{2} recorded_LS{3} recorded_LS{4} recorded_LS{5}];

%find unique features
unique_features = unique(feature_combination);

%count how many times each unique feature is selected over five folds in
%the inner loop
counts = [];
for i = 1:length(unique_features)
counts = [counts sum(feature_combination==unique_features(i))];
end

%calculate the average log length scale of each unique features over five folds in
%the inner loop
LSvalues = [];
for i = 1:length(unique_features)
LSvalues = [LSvalues sum(LS_combination(feature_combination==unique_features(i)))/sum(feature_combination==unique_features(i))];
end

%sort unique features and their corresponding average log length scales based on 
%their number of hits
[sorted_counts index_counts] = sort(counts,'descend');
unique_features_sorted = unique_features(index_counts);
LSvalues_sorted = LSvalues(index_counts);

%resort unique features with the same number of hits based on their average
%log length scales
resorted_featurecounts = [];
resorted_LSvalues = [];
for i = 1:5
    subset = find(sorted_counts==6-i);
    sub_table = [unique_features_sorted(subset)' sorted_counts(subset)' LSvalues_sorted(subset)'];

    [LS_subset index_subset] = sort(LSvalues_sorted(subset));
    sub_table = [sub_table(index_subset,1) sub_table(index_subset,2)];
    
    resorted_LSvalues = [resorted_LSvalues; LS_subset'];
    resorted_featurecounts = [resorted_featurecounts; sub_table];
    
end

%find the resorted unique features
sorted_unique_features = resorted_featurecounts(:,1); 

%carry out forward selection to calclulate LOO errors again using the
%selected features
final_set = train_valid_set(:,sorted_unique_features);
labels = train_valid_label;
LOOerr = forwardFS_LOO(final_set,labels);

%find the minimal subset of features achieving the lowest LOO error
min_loc = find(LOOerr==min(LOOerr))-1;
final_features = sorted_unique_features(1:min_loc);

%reduce the train-valid data
final_trainvalid = train_valid_set(:,final_features);
final_trainvalid_label = train_valid_label;

%reduce the test data
final_test = test_set(:,final_features);
final_test_label = test_label;

%find the rank of the selected features in t-test
rank_t = [];
for i = 1:length(final_features)
r = find(IDX_t==final_features(i));
rank_t = [rank_t r];
end

%calculate LOO errors using top-ranked 50 features selected by t-test 
final_set2 = train_valid_set(:,IDX_t(1:50));
labels2 = train_valid_label;
LOOerr2 = forwardFS_LOO(final_set2,labels2);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set the log initial values of two hyperparameters (a length scale and a signal
%variance)
%set the length scale to be the average of pairwise Euclidean distances 
avgD = mean(mean(pdist2(final_trainvalid,final_trainvalid)));
loghyper = [log(avgD); 0];

%obtain the predictive probabilities using initial hyperparameter values    
p = binaryLaplaceGP(loghyper, 'covSEiso', 'cumGauss',final_trainvalid,final_trainvalid_label, final_test);
  
%claculate prediction accuracy and information in bits  
accuracy1 = 1-sum((p>0.5)~=(final_test_label>0))/length(final_test_label);
imf1 = mean((final_test_label==1).*log2(p)+(final_test_label==-1).*log2(1-p))+1;
  
%maximize the log marginal likelihood w.r.t hyperparameters by carrying out
%conjugate gradient optimization (-50 tells minimize to evaluate the function at most 50 times)  
[loghyper' binaryLaplaceGP(loghyper, 'covSEiso', 'cumGauss',final_trainvalid, final_trainvalid_label)]
[newloghyper logmarglik] = minimize(loghyper, 'binaryLaplaceGP', -50, 'covSEiso', 'cumGauss', final_trainvalid, final_trainvalid_label);
[newloghyper' logmarglik(end)]
  
%obtain the predictive probabilities using optimized hyperparameter values   
pp = binaryLaplaceGP(newloghyper, 'covSEiso', 'cumGauss', final_trainvalid, final_trainvalid_label, final_test);
  
%claculate new prediction accuracy and information in bits    
accuracy2 = 1-sum((pp>0.5)~=(final_test_label>0))/length(final_test_label);
imf2 = mean((final_test_label==1).*log2(pp)+(final_test_label==-1).*log2(1-pp))+1;
 
%draw a plot of predictive probabilities calculated by using both initial
%and optimized hyperparameter values 
figure;
h = plot(p(1:16),'.','MarkerSize',14,'MarkerEdgeColor',[0,0.7,0.8],'MarkerFaceColor',[0,0.7,0.8])
hold on  
plot(17:32,p(17:32),'^','MarkerSize',4,'MarkerEdgeColor',[0,0.7,0.8],'MarkerFaceColor',[0,0.7,0.8])
plot(0:35,repmat(0.5,1,36),'r') 
plot(pp(1:16),'.','MarkerSize',14,'MarkerEdgeColor',[0,0.3,0.4],'MarkerFaceColor',[0,0.3,0.4])
plot(17:32,pp(17:32),'^','MarkerSize',4,'MarkerEdgeColor',[0,0.3,0.4],'MarkerFaceColor',[0,0.3,0.4]) 
legend('AD (initial hyperparameter values)','Control (initial hyperparameter values)','Thresholding the predictive probabilities at 0.5','AD (optimized hyperparameter values)','Control (optimized hyperparameter values)')  
title('Test set predictions')
xlabel('Test case label') 
ylabel('Predictive probability') 
hold off 

%save the figure
saveas(h,fullfile('C:\Work Files\project\data\',[sprintf('FIG%d',ii) '.jpeg']));

%record both the prediction accuracies and the information in bits
result = [accuracy1; imf1; accuracy2; imf2];
output = [output result];

%record both the LOO errors calculated by using 50 features sorted by the
%number of hits and 50 top-ranked features in ttest
LOOE{ii} = LOOerr;
LOOE2{ii} = LOOerr2;

%record average log length scales
Avg_log_LS{ii} = resorted_LSvalues;

%record the selected features and their number of hits
Features_Hits{ii} = resorted_featurecounts;

%record the rank of the selected features in ttest
Rank_in_ttest{ii} = rank_t;

end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate the mean and standard deviation of the initial accuracies
initial_accuracies = output(1,:);
mean(initial_accuracies)
std(initial_accuracies,[],2)

%calculate the mean and standard deviation of the optimized accuracies
optimized_accuracies = output(3,:);
mean(optimized_accuracies)
std(optimized_accuracies,[],2)






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%draw a plot of LOO errors (using 50 features sorted by the
%number of hits) obtained in the first fold in the outer loop
LOOerr = LOOE{1};
plot(1:50,LOOerr(2:51),'o')
ylim([0 50])
title('Leave-One-Out error using the top-ranked features')
xlabel('Number of features in modelling') 
ylabel('Leave-One-Out error')
hold on
plot(5,14,'.','MarkerSize',20)

%draw a plot of LOO errors (using 50 top-ranked features in ttest)obtained 
%in the first fold in the outer loop
LOOerr2 = LOOE2{1};
plot(1:50,LOOerr2(2:51),'o')
ylim([0 50])
title('Leave-One-Out error using the top-ranked features')
xlabel('Number of features in modelling') 
ylabel('Leave-One-Out error')
hold on
plot(11,19,'.','MarkerSize',20)






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find the names of the selected genes and some random genes
%'a' and 'index_delete' are obtained by running the code given in the file 
%called 'preprocess_visualization'
allgenes = a(1:48228);
genes = allgenes(~index_delete(1:48228));
genelabels = genes(final_features);
genelabels2 = genes(final_features+990);

%find sample IDs labelled by classes (ADs or Controls)
sampleID = cell(160,1);
for i = 1:82
sampleID{i} = ['AD' b{i}];
end
for i = 83:160
sampleID{i} = ['Ctrl' b{i}];
end

%draw clustering trees and heat map using features selected in the first
%fold in the outer loop
heatmap_data = transpose(test_set(:,final_features));
cgo = clustergram(heatmap_data,'Standardize','none','RowPDist', {'correlation'},'ColumnPDist',{'correlation'})
set(cgo,'RowLabels',cellstr(genelabels),'ColumnLabels',cellstr(sampleID(Indices2==1))')

%draw clustering trees and heat map using random features 
%the number of these random features is the same as that of the selected
%features
heatmap_data = transpose(test_set(:,final_features+990));
cgo = clustergram(heatmap_data,'Standardize','none','RowPDist', {'correlation'},'ColumnPDist',{'correlation'})
set(cgo,'RowLabels',cellstr(genelabels2),'ColumnLabels',cellstr(sampleID(Indices2==1))')






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%draw correlation heat map using all features
figure;
%colormap('hot')
imagesc(corrcoef(tb_aftctl'))
colorbar

%draw correlation heat map using 100 top-ranked features in ttest
figure;
%colormap('hot')
imagesc(corrcoef(tb_aftctl(:,IDX_t)'))
colorbar

%draw correlation heat map using features selected by using the proposed
%resampling procedure
figure;
%colormap('hot')
imagesc(corrcoef(tb_aftctl(:,final_features)'))
colorbar

