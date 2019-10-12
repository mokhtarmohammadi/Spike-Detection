 FPR=mean(FPR1'/100);
 TPR=mean(Hitrate1'/100);
 k = length(FPR);
 A = zeros(k,1);

if length(TPR) ~= k
    disp('Length of vectors FPR and TPR must be equal.');
    return;
else
    figure;
    spec=FPR;
    sens=TPR;
     cspec = 1-spec;
cspec = cspec(end:-1:1);
sens = sens(end:-1:1);
plot(cspec,sens,'k-')
     set(gca,'FontSize',25);
AUC = sum(0.5*(sens(2:end)+sens(1:end-1)).*(cspec(2:end) - cspec(1:end-1)));
fprintf('\nAUC: %g',AUC);
  
end
%  for i = 2 : k
%     A(i) = (FPR(i - 1) - FPR(i)).*(1 - TPR(i));
% end
% area1 = 1 - sum(A)

%%%%%%%%%%
 FPR=mean(FPR2'/100);
 TPR=mean(Hitrate2'/100);
 k = length(FPR);
 A = zeros(k,1);

if length(TPR) ~= k
    disp('Length of vectors FPR and TPR must be equal.');
    return;
else
  hold on;
    spec=FPR;sens=TPR;
     cspec = 1-spec;
cspec = cspec(end:-1:1);
sens = sens(end:-1:1);
plot(cspec,sens,'k-')
     set(gca,'FontSize',25);
AUC = sum(0.5*(sens(2:end)+sens(1:end-1)).*(cspec(2:end) - cspec(1:end-1)));
fprintf('\nAUC: %g',AUC);
  
end

% % for i = 2 : k
% %     A(i) = (FPR(i - 1) - FPR(i)).*(1 - TPR(i));
% % end
% area2 = 1 - sum(A) 

%%%%%%%%%%%%%%
 FPR=mean(FPR3'/100);
 TPR=mean(Hitrate3'/100);
 k = length(FPR);
 A = zeros(k,1);

if length(TPR) ~= k
    disp('Length of vectors FPR and TPR must be equal.');
    return;
else
    hold on;
    spec=FPR;sens=TPR;
     cspec = 1-spec;
cspec = cspec(end:-1:1);
sens = sens(end:-1:1);
plot(cspec,sens,'k-')
     set(gca,'FontSize',25);
AUC = sum(0.5*(sens(2:end)+sens(1:end-1)).*(cspec(2:end) - cspec(1:end-1)));
fprintf('\nAUC: %g',AUC);
  
%     axis([0 1 0 1]);
end
% for i = 2 : k
%     A(i) = (FPR(i - 1) - FPR(i)).*(1 - TPR(i));
% end
% area3 = 1 - sum(A) 
 title('ROC Curve');
 xlabel('False Positive Rate(in dB)');
 ylabel('Hit Rate');
hold off;