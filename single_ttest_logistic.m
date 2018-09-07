%replace the original class lablels (1 and 0) by new ones (1 and -1)
label_set = tb_aftctl(:,47202);
label_set = sign(label_set-0.5);

%generate the indices for five-fold cross-validation 
Indices = crossvalind('Kfold',label_set, 5);

sd_all = [];
m_all = [];
 
for j = [47201 20000 10000 5000 1000 500 200 75 25] 

output = [];
  
for i = 1:5
    
%four folds are used as training data
train_set = tb_aftctl(Indices~=i,1:47201);
train_label = label_set(Indices~=i);

%one fold left out is regarded as test dataset
test_set = tb_aftctl(Indices==i,1:47201);
test_label = label_set(Indices==i);

%perform t-test within the training data (assume that the variances of different groups are unequal)
[h,pvalue,ci,stat] = ttest2(train_set(train_label==1,:),train_set(train_label==-1,:),'Vartype','unequal');
 
%find indices that sorts the p-values in ascending order 
[sorted_pvalue,index] = sort(pvalue);
 
%select the first j features with smaller p-values
IDX = index(1:j); 

%reduce both training and test datasets 
train_set = train_set(:,IDX);
test_set = test_set(:,IDX); 

%perform Gaussian process classification using isotropic squared
%exponential covariance function and logistic likelihood function

%set the log initial values of two hyperparameters (a length scale and a signal
%variance)
%set the length scale to be the average of pairwise Euclidean distances 
avgD = mean(mean(pdist2(train_set,train_set)));

%set the signal variance to be 1
loghyper = [log(avgD); 0];
  
%obtain the predictive probabilities using initial hyperparameter values    
p = binaryLaplaceGP(loghyper, 'covSEiso', 'logistic', train_set, train_label, test_set);
  
%claculate prediction accuracy and information in bits 
accuracy1 = 1-sum((p>0.5)~=(test_label>0))/length(test_label);
imf1 = mean((test_label==1).*log2(p)+(test_label==-1).*log2(1-p))+1;
  
%maximize the log marginal likelihood w.r.t hyperparameters by carrying out
%conjugate gradient optimization (-50 tells minimize to evaluate the function at most 50 times)
[loghyper' binaryLaplaceGP(loghyper, 'covSEiso', 'logistic', train_set, train_label)]
[newloghyper logmarglik] = minimize(loghyper, 'binaryLaplaceGP', -50, 'covSEiso', 'logistic', train_set, train_label);
[newloghyper' logmarglik(end)]
  
%obtain the predictive probabilities using optimized hyperparameter values 
pp = binaryLaplaceGP(newloghyper, 'covSEiso', 'logistic', train_set, train_label, test_set);
  
%claculate new prediction accuracy and information in bits   
accuracy2 = 1-sum((pp>0.5)~=(test_label>0))/length(test_label);
imf2 = mean((test_label==1).*log2(pp)+(test_label==-1).*log2(1-pp))+1;
  
%record the results for a particular number of features 
result = [accuracy1; imf1; accuracy2; imf2];
output = [output result];
  
end

%find the mean and the standard deviation of the prediction accuracies (as well as those of the 
%information in bits) over five folds
m=mean(output,2);
sd=std(output,[],2);

%record the mean and the standard deviation
m_all = [m_all m];
sd_all = [sd_all sd];

end
  
%view the results obtained by selecting different number of features
m_all 
sd_all 
  





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create bar plots of prediction accuracies 
figure
acc = [m_all(1,:)' m_all(3,:)'];
h = bar(acc)
xticklabels({'47201', '20000', '10000', '5000', '1000', '500', '200', '75', '25'})
ylim([0 1])
colormap(summer(n));
grid on
hold on;
  
%add a horizontal line reflecting majority class prediction
base = 82/160;
hline = refline([0 base]);
hline.Color = 'r';
hline.LineWidth = 1.5;
 
%add a legend
leg1 = legend({'Prediction using initial hyperperameter values','Prediction using optimized hyperperameter values','Majority class prediction'},'FontSize',12)
  
%add error bars reflecting the standard deviations of the prediction
%accuracies
x = m_all(1,:);
y = sd_all(1,:);
errorbar([1:9]-0.14,x,y,'.k','LineWidth',1.5,'HandleVisibility','off')  
x = m_all(3,:);
y = sd_all(3,:);
errorbar([1:9]+0.14,x,y,'.k','LineWidth',1.5,'HandleVisibility','off')
 
%add titles and axis labels
title('Prediction Accuracy (T-test)')
xlabel('Number of features') 
ylabel('Accuracy')
 
hold off
  
%create bar plots of information in bits
figure
acc = [m_all(2,:)' m_all(4,:)'];
h = bar(acc)
xticklabels({'47201', '20000', '10000', '5000', '1000', '500', '200', '75', '25'})
ylim([0 0.5])
colormap(summer(n));
grid on
hold on;
 
%add a legend
legend({'Prediction using initial hyperperameter values','Prediction using the optimized hyperperameter values'},'FontSize',12)
  
%add error bars reflecting the standard deviations of the information in bits 
x = m_all(2,:);
y = sd_all(2,:);
errorbar([1:9]-0.14,x,y,'.k','LineWidth',1.5,'HandleVisibility','off') 
x = m_all(4,:);
y = sd_all(4,:);
errorbar([1:9]+0.14,x,y,'.k','LineWidth',1.5,'HandleVisibility','off')

%add titles and axis labels
title('Information in bits about the targets (T-test)')
xlabel('Number of features') 
ylabel('Information')
  
hold off