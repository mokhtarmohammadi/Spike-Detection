function [area] = rocplot(FPR,TPR)

% Usage [area] = rocplot(FPR,TPR)
% It plots the ROC curve as TPR vs. FPR
% where TPR is the true positive rate and
% FPR is the false positive rate. TPR and
% FPR must be vectors of equal size. For
% detail see Chaovalitwongse, Iasemidis et
% al. "Performance of a seizure warning
% algorithm based on the dynamics of
% intracranial EEG," Epilepsy Res. 64 : 93
% - 113, 2005.
% The program also calculates an 
% approximate area under the ROC curve. The
% higher the number of data points the more
% accurate will be the claculated area.
% Copy right: Kaushik Majumdar, Oct 2009.

k = length(FPR);
A = zeros(k,1);

if length(TPR) ~= k
    disp('Length of vectors FPR and TPR must be equal.');
    return;
else
    plot(FPR,TPR);
    line([0 FPR(k)], [0 TPR(k)]);
    line([1 FPR(1)], [1 TPR(1)]);
    axis([0 1 0 1]);
    title('ROC Curve');
    xlabel('FPR');
    ylabel('TPR');
end

for i = 2 : k
    A(i) = (FPR(i - 1) - FPR(i)).*(1 - TPR(i));
end

area = 1 - sum(A);