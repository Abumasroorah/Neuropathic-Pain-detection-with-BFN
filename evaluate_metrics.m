function [acc, prec, rec, f1] = evaluate_metrics(y_true, y_pred)
    % Confusion matrix
    cm = confusionmat(y_true, y_pred);
    num_classes = size(cm, 1);

    % True Positives, False Positives, False Negatives
    TP = diag(cm);
    FP = sum(cm, 1)' - TP;
    FN = sum(cm, 2) - TP;

    % Per-class metrics
    precision = TP ./ (TP + FP + eps);
    recall = TP ./ (TP + FN + eps);
    f1_class = 2 * (precision .* recall) ./ (precision + recall + eps);

    % Weighted average
    support = sum(cm, 2);  % Samples per class
    total = sum(support);
    acc = sum(TP) / total;
    prec = sum(support .* precision) / total;
    rec = sum(support .* recall) / total;
    f1 = sum(support .* f1_class) / total;
end
