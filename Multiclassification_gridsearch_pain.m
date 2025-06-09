clear;
close all;
clc;
load New_featureEC2;
% Load the data
% data = load('FTEO_NCHPMPLP.mat');
% X = data.ft_Gamma(:, [1:3 5]);  % Select 4 features
% X = ft_Gamma(:, [1:3 5]);
X = Feat_matrix_delt(:, [1:3 5:10]);

% Compute wavelet coefficients (e.g., 4-40 Hz)
X_aug = [X, log(X + eps), sqrt(X), X.^2];  % Expanded features
[n_samples, n_features] = size(X);
X_pair = [];
for i = 1:n_features
    for j = i+1:n_features
        X_pair = [X_pair, X(:, i) .* X(:, j)];  % element-wise product
    end
end
X_aug = [X_aug, X_pair];

% y = double(data.Lab(:));        % Ensure y is a column vector
y = Label_vector;        % Ensure y is a column vector

% Shuffle the data
rng(42);  % For reproducibility
idx = randperm(size(X,1));
X = X(idx, :);
y = y(idx);

% Stratified split: 80% train, 10% val, 10% test
cv = cvpartition(y, 'HoldOut', 0.2);
X_train_full = X(training(cv), :);
y_train_full = y(training(cv));
X_temp = X(test(cv), :);
y_temp = y(test(cv));

cv2 = cvpartition(y_temp, 'HoldOut', 0.5);
X_val = X_temp(test(cv2), :);
y_val = y_temp(test(cv2));
X_test = X_temp(training(cv2), :);
y_test = y_temp(training(cv2));

% Define a compact manual grid
hyperparams = {
    struct('BoxConstraint', 0.1, 'KernelFunction', 'linear'), ...
    struct('BoxConstraint', 1,   'KernelFunction', 'linear'), ...
    struct('BoxConstraint', 0.1, 'KernelFunction', 'gaussian', 'KernelScale', 1), ...
    struct('BoxConstraint', 1,   'KernelFunction', 'gaussian', 'KernelScale', 1), ...
    struct('BoxConstraint', 1,   'KernelFunction', 'gaussian', 'KernelScale', 0.5), ...
    struct('BoxConstraint', 0.1, 'KernelFunction', 'polynomial', 'PolynomialOrder', 2), ...
    struct('BoxConstraint', 1,   'KernelFunction', 'polynomial', 'PolynomialOrder', 2), ...
    struct('BoxConstraint', 0.1, 'KernelFunction', 'polynomial', 'PolynomialOrder', 3), ...
    struct('BoxConstraint', 1,   'KernelFunction', 'polynomial', 'PolynomialOrder', 3)
};


bestCVLoss = Inf;
bestTemplate = [];
bestDesc = '';

% Grid search with 5-fold cross-validation
for i = 1:length(hyperparams)
    hp = hyperparams{i};
    
    % Default values
    ks = 'auto';
    po = [];
    
    if isfield(hp, 'KernelScale')
        ks = hp.KernelScale;
    end
    if isfield(hp, 'PolynomialOrder')
        po = hp.PolynomialOrder;
    end
    
    % Create SVM template
    t = templateSVM( ...
        'KernelFunction', hp.KernelFunction, ...
        'BoxConstraint', hp.BoxConstraint, ...
        'KernelScale', ks, ...
        'PolynomialOrder', po ...
    );

    % 5-fold cross-validation with ECOC
    Mdl = fitcecoc(X_train_full, y_train_full, ...
        'Learners', t, ...
        'Coding', 'onevsall', ...
        'KFold', 5);

    cvLoss = kfoldLoss(Mdl);
    fprintf('Model %d: Kernel=%s, C=%.2f, Loss=%.4f\n', i, hp.KernelFunction, hp.BoxConstraint, cvLoss);
    
    % Save best model
    if cvLoss < bestCVLoss
        bestCVLoss = cvLoss;
        bestTemplate = t;
        bestDesc = sprintf('Kernel=%s, C=%.2f', hp.KernelFunction, hp.BoxConstraint);
    end
end

% Train final model with best hyperparameters
finalModel = fitcecoc(X_train_full, y_train_full, ...
    'Learners', bestTemplate, ...
    'Coding', 'onevsall');

% Predict on test and validation sets
y_pred_test = predict(finalModel, X_test);
y_pred_val = predict(finalModel, X_val);


% Evaluate performance
[accuracy_test, precision_test, recall_test, f1_test] = evaluate_metrics(y_test, y_pred_test);
[accuracy_val, precision_val, recall_val, f1_val] = evaluate_metrics(y_val, y_pred_val);

% Display results
fprintf('\nBest Model: %s\n', bestDesc);
fprintf('Testing Set Performance:\n');
fprintf('Accuracy: %.4f\n', accuracy_test);
fprintf('Precision: %.4f\n', precision_test);
fprintf('Recall: %.4f\n', recall_test);
fprintf('F1-Score: %.4f\n\n', f1_test);

fprintf('Validation Set Performance:\n');
fprintf('Accuracy: %.4f\n', accuracy_val);
fprintf('Precision: %.4f\n', precision_val);
fprintf('Recall: %.4f\n', recall_val);
fprintf('F1-Score: %.4f\n', f1_val);
