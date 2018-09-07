function f = gp_LOO(train_set,train_label,test_set,test_label,cov,likeli)
  %this function is only used in the 'forwardFS_LOO' function to simplify
  %the code
  avgD = mean(mean(pdist2(train_set,train_set)));
  loghyper = [log(avgD); 0];

  p = binaryLaplaceGP(loghyper, cov, likeli, train_set, train_label, test_set);
  
  [loghyper' binaryLaplaceGP(loghyper, cov, likeli, train_set, train_label)];
  [newloghyper logmarglik] = minimize(loghyper, 'binaryLaplaceGP', -50,cov,likeli, train_set, train_label);
  [newloghyper' logmarglik(end)];
  
  pp = binaryLaplaceGP(newloghyper, cov, likeli, train_set, train_label, test_set);
    
  f = sum((pp>0.5)~=(test_label>0));
  
end