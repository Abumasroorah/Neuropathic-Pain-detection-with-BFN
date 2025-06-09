clear;
close all;
clc;

load featureEC;

% Feature selection and augmentation
% X = Feat_matrix_delt(:, [1:3 5:10]);
X= X_HP_vs_LP(:, [1:3 5:10]);
X_aug = [X, log(X + eps), sqrt(X), X.^2];  % Augmented features
% y = Label_vector;
y=double(y_HP_vs_LP);

% Shuffle data
rng(42);
idx = randperm(size(X_aug,1));
X = X_aug(idx, :);
y = y(idx);

% Define a compact manual grid of hyperparameters
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

bestF1 = -inf;
bestParams = struct();
bestDesc = '';

% Grid search with 5-fold cross-validation
for i = 1:length(hyperparams)
    hp = hyperparams{i};

    % Set default optional parameters
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

    % 5-fold CV using ECOC model
    Mdl = fitcecoc(X, y, 'Learners', t, 'Coding', 'onevsall', 'KFold', 5);

    % Collect predictions and true labels from all folds
    y_cv_true = [];
    y_cv_pred = [];

    for k = 1:Mdl.KFold
        trainedModel = Mdl.Trained{k};
        testIdx = Mdl.Partition.test(k);
        X_cv = X(testIdx, :);
        y_cv = y(testIdx);

        y_pred_k = predict(trainedModel, X_cv);
        y_cv_true = [y_cv_true; y_cv];
        y_cv_pred = [y_cv_pred; y_pred_k];
    end

    % Evaluate metrics
    [acc_cv, prec_cv, rec_cv, f1_cv] = evaluate_metrics(y_cv_true, y_cv_pred);

    % Display fold performance
    fprintf('Model %d: Kernel=%s, C=%.2f | CV Acc: %.4f | Precision: %.4f | Recall: %.4f | F1: %.4f\n', ...
        i, hp.KernelFunction, hp.BoxConstraint, acc_cv, prec_cv, rec_cv, f1_cv);

    % Save best
    if f1_cv > bestF1
        bestF1 = f1_cv;
        bestParams = hp;
        bestDesc = sprintf('Kernel=%s, C=%.2f', hp.KernelFunction, hp.BoxConstraint);
    end
end

% Retrain final model on all data using best hyperparameters
ks = 'auto'; po = [];
if isfield(bestParams, 'KernelScale'), ks = bestParams.KernelScale; end
if isfield(bestParams, 'PolynomialOrder'), po = bestParams.PolynomialOrder; end

finalModel = fitcsvm(X, y, ...
    'KernelFunction', bestParams.KernelFunction, ...
    'BoxConstraint', bestParams.BoxConstraint, ...
    'KernelScale', ks, ...
    'PolynomialOrder', po);

% Final model info
fprintf('\nBest Model Selected:\n%s | Final F1 (CV): %.4f\n', bestDesc, bestF1);
