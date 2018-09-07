%find gene expression data where the effects from sex, age and RNA quality
%have been removed 
sett = tb_aftctl(:,1:47201);

%find the sample labels (ADs are labelled by +1 while Contrals are labelled
%by -1)
label_set = tb_aftctl(:,47202);
label_set = sign(label_set-0.5);
label = label_set;

%perform ttest based on the whole data (the whole data can be used here
%since we do not aim to train any classifier but just check feature importances) 
[h,pvalue,ci,stat] = ttest2(sett(label==1,1:47201),sett(label==-1,:),'Vartype','unequal');

%draw the cumulative density function of pvalues
ecdf(pvalue);
hold on

%find a point indicating the percentages of features with p-values less
%than 0.01 in the CDF plot 
scatter(0.01,sum(pvalue<=0.01)/47201,'filled')
x=0:0.01:1;
y=sum(pvalue<=0.01)/47201;
plot(x,y*ones(size(x)),'-.r','HandleVisibility','off')
x=0.01;
y=0:0.01:1;
plot(x*ones(size(y)),y,'-.r','HandleVisibility','off')

%find a point indicating the percentages of features with p-values less
%than 0.05 in the CDF plot 
scatter(0.05,sum(pvalue<=0.05)/47201,'filled')
x=0:0.01:1;
y=sum(pvalue<=0.05)/47201;
plot(x,y*ones(size(x)),'-.m','HandleVisibility','off')
x=0.05;
y=0:0.01:1;
plot(x*ones(size(y)),y,'-.m','HandleVisibility','off')

%add a legend
legend('CDF Curve','Point (0.01, 0.18)','Point (0.05, 0.27)')

%%add titles and axis labels
title('Cumulative density function plot')
xlabel('P value');
ylabel('CDF value');

hold off