function E = forwardFS_LOO(data_set,data_set_label)
% this function is used to calculate LOO errors by carrying out forward
% selection (a feature is added to the modelling at a time)
E = size(data_set,1);
    
for j = 1:size(data_set,2)
    
X = data_set(:,1:j);
y = data_set_label;

err = 0;

%perform Leave-One-Out cross-validation given dataset with j features
for i = 1:size(data_set,1)
    
    index_LOO = 1:size(data_set,1);
    
    X_train = X(index_LOO~=i,:);
    y_train = y(index_LOO~=i);
    X_valid = X(index_LOO==i,:);
    y_valid = y(index_LOO==i);
    
    %count the number of errors
    err = err + gp_LOO(X_train,y_train,X_valid,y_valid,'covSEiso','cumGauss');
    
    disp([i j])
end

E = [E err];

end



