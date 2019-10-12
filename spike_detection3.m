close all;
clear all;
addpath('D:\tfsa_5-5\windows\win64_bin');
x=50*cos(2*pi*0.05*(0:255)+0*2*pi*0.0000007*(0:255).^3)+100*cos(2*pi*0.4*(0:255)-2*pi*0.0000007*(0:255).^3);
 x_s=zeros(1,256);
 %x_s(40:80:end)=4;
 x_s(40)=1;
 x_s(41)=0;
 x_s(42)=0.0;
 x_s(39)=0;
 x_s(38)=0.0;
 
 
 x_s(120-2:120+2)=x_s(40-2:40+2);
 x_s(200-2:200+2)=x_s(40-2:40+2);
 
 x=x+100*x_s;
 x=x.';
 fs=1;
 mintime_m=0.003;
 SNR=[-15 -10 -5 0 5 10 15];
 for nn=7:7
    for mmm=1:10
    x=awgn(x,SNR(nn),'measured');
    end
     actualpeaks=[40 120 200];
    [adtfd1,orient]=HTFD_new2_spike(x,3,5,100);
    
    Mask=zeros(size(orient));
     orient=(orient-1)*45;
    Mask( and(orient<110,orient>70))=1;
    
    Mask(:,1:15)=0;
    Mask(:,end-15:end)=0;
    
    adtfd=Mask.*adtfd1;
    %tfsapl(x,adtfd)
    %Is=quadtfd(x,length(x)-1,1,'specx',15,'hamm',256);
    %figure; tfsapl(x,adtfd)
    %I2=II;
    %II(I2<0.2*max(max(II)))
    MASK=adtfd/max(adtfd(:));
    MASK(MASK>0.1)=1;
    MASK(MASK<0.1)=0;
    
%     se = strel('line',10,90);
%     
%     MASK=imerode(MASK,se);
%     se = strel('line',10,90);
%     
%     MASK=imdilate(MASK,se);
%     
%     
    
    
   % mm=sum(adtfd.*MASK);
    mm=sum(adtfd);
    mm=mm/max(mm);
    mm(mm<0.2)=0;
     %Estimate the location of the spikes in Seconds
    [p1, kest1]=findpeaks(mm,'MINPEAKHEIGHT',0.3);
%     if ~isempty(spikes); c_spike = zeros(length(dataset), 1); c_spike (spikes)= 1; end;
%     kest1=find(c_spike(1:end)>0);
    id=find(diff(kest1)< (mintime_m*fs));
    if ~isempty(id)
        temp=(1:id(1)); for i=2:length(id) temp=[temp 2+id(i-1):id(i)]; end
        temp=[temp 2+id(end): length(kest1)]; kest=kest1(temp);
    else kest=kest1;
    end
    pos_len1 = length(kest);
    T_pos1 = zeros(1, pos_len1);
    for i=1:pos_len1
        T_pos1(i) = kest(i) /fs;
    end
    % find number of missing spikes
    [TP1 FN1] = matchspikes_sd(actualpeaks, T_pos1,'MinTime', mintime_m, 'SampleRate', (fs)) ;
    %find inserted spikes
% FP1 = findinsertions_sd2(actualpeaks, T_pos1,'MinTime', mintime_m) ;
    FP1=length(kest)-TP1;
    TN1=length(x)- length(actualpeaks)-FN1;
    %%%%%%%%%%%%%%%%
    Hitrate1(nn,mmm) = (TP1/(TP1+FN1))*100;
    Precision1(nn,mmm)=( TP1/(TP1+FP1))*100;
    FPR1(nn,mmm)= (FP1/(TN1+FP1))*100;
    ATP1(nn,mmm)=(TP1/length(actualpeaks))*100;
    AFP1(nn,mmm)=(FP1/length(actualpeaks))*100;
    %%%%%%%%%%%%%%%%%
    figure;plot(kest,1,'*r');
    title('ADTFD','FontSize',20);
    hold on;
    plot(mm);set(gca,'FontSize',20);
     cc=round(actualpeaks*fs);
    
%     cc=cc(cc<length(x)-10);
    
    hold on;plot(cc,1.1,'Ok')
    hold off;
    clear temp kest kest1;
    %%%%%%%%%%%%%%%%%%%%%%%  SNEO  %%%%%%%%%%%
    yy=x(2:end-1).^2-x(3:end).*x(1:end-2);
    %yy=filter(gausswin(20,1),1,yy);
    %yy=[0;0; yy]';
    yy=yy/max(abs(yy));
    yy1=find(yy>0.2*max(abs(yy)) );
    [p1, kest1]=findpeaks(yy,'MINPEAKHEIGHT',0.3) ;
    
    id=find(diff(kest1)< (mintime_m*fs));
    if ~isempty(id)
        temp=(1:id(1)); for i=2:length(id) temp=[temp 2+id(i-1):id(i)]; end
        temp=[temp 2+id(end): length(kest1)]; kest=kest1(temp);
    else kest=kest1;
    end
%     kest=kest1;
    pos_len1 = length(kest);
    T_pos1 = zeros(1, pos_len1);
    for i=1:pos_len1
        T_pos1(i) = kest(i) /fs;
    end
    % find number of missing spikes
    [TP2 FN2] = matchspikes_sd(actualpeaks, T_pos1,'MinTime', mintime_m, 'SampleRate', (fs)) ;
    %find inserted spikes
% FP2 = findinsertions_sd2(actualpeaks, T_pos1,'MinTime', mintime_m) ;
    FP2=length(kest)-TP2;
    TN2=length(x)- length(actualpeaks)-FN2;
    %%%%%%%%%%%%%%%%
    Hitrate2(nn,mmm) = (TP2/(TP2+FN2))*100;
    Precision2(nn,mmm)=( TP2/(TP2+FP2))*100;
    FPR2(nn,mmm)= (FP2/(TN2+FP2))*100;
    ATP2(nn,mmm)=(TP2/length(actualpeaks))*100;
    AFP2(nn,mmm)=(FP2/length(actualpeaks))*100;
    %%%%%%%%%%%%%%%%%
%     figure; plot(yy);set(gca,'FontSize',20);
%     hold on;
%     plot(kest,1,'Xr');
%     %
%     title('SNEO');
%     hold on;plot(cc,1.1,'Og');
%     hold off;
    clear temp kest kest1;
    %%%%%%%%%%%%%%%%%    COB    %%%%%%%
    peaktrain=spbycob(x, 64);
    %=========   thresholding ===========
    peaktrain=peaktrain/max(abs(peaktrain));
    [p3, kest1]=findpeaks(peaktrain,'MINPEAKHEIGHT',0.3) ;
    id=find(diff(kest1)< (mintime_m*fs));
    if ~isempty(id)
        temp=(1:id(1)); for i=2:length(id) temp=[temp 2+id(i-1):id(i)]; end
        temp=[temp 2+id(end): length(kest1)]; kest=kest1(temp);
    else kest=kest1;
    end
%     kest=kest1;
    pos_len1 = length(kest);
    T_pos1 = zeros(1, pos_len1);
    for i=1:pos_len1
        T_pos1(i) = kest(i) /fs;
    end
    % find number of missing spikes
    [TP3 FN3] = matchspikes_sd(actualpeaks, T_pos1,'MinTime', mintime_m, 'SampleRate', (fs)) ;
    %find inserted spikes
% FP3 = findinsertions_sd2(actualpeaks, T_pos1,'MinTime', mintime_m) ;
    FP3=length(kest)-TP3;
    TN3=length(x)- length(actualpeaks)-FN3;
    %%%%%%%%%%%%%%%%
    Hitrate3(nn,mmm) = (TP3/(TP3+FN3))*100;
    Precision3(nn,mmm)=( TP3/(TP3+FP3))*100;
    FPR3(nn,mmm)= (FP3/(TN3+FP3))*100;
    ATP3(nn,mmm)=(TP3/length(actualpeaks))*100;
    AFP3(nn,mmm)=(FP3/length(actualpeaks))*100;
    %%%%%%%%%%%%%%%%%
%     figure;plot(peaktrain);set(gca,'FontSize',20);
%     title('COB');
%     hold on;
%     plot(kest,1,'*r');
%     hold on;plot(cc,1.1,'Og')
%     hold off;
    clear detected temp kest1 kest2 kest3;
    end

figure;
plot(SNR,(mean(Hitrate1')),'k-o','MarkerSize',10,'LineWidth',3)
 set(gca, 'XDir','reverse','FontSize',25);
 hold on;
 plot(SNR,(mean(Hitrate2')),'k--s','MarkerSize',10,'LineWidth',3)
 hold on;
 plot(SNR,(mean(Hitrate3')),'k:>','MarkerSize',10,'LineWidth',3)
 hold off;
 xlabel('SNR (in dB)');
 ylabel('Hit Rate (in %)');
 legend('ADTFD','PSNEO','COB');
 %%%%%%%%%%%%%%%%%
 figure;
 plot(SNR,(mean(Precision1')),'k-o','MarkerSize',10,'LineWidth',3)
  set(gca, 'XDir','reverse','FontSize',25);
 hold on;
 plot(SNR,(mean(Precision2')),'k--s','MarkerSize',10,'LineWidth',3)
 hold on;
 plot(SNR,(mean(Precision3')),'k:>','MarkerSize',10,'LineWidth',3)
 hold off;
 xlabel('SNR (in dB)');
 ylabel('Precision (in %)');
  legend('ADTFD','PSNEO','COB');
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 figure;
plot(SNR,(mean(ATP1')),'k-o','MarkerSize',10,'LineWidth',3)
 set(gca, 'XDir','reverse','FontSize',25);
 hold on;
 plot(SNR,(mean(ATP2')),'k--s','MarkerSize',10,'LineWidth',3)
 hold on;
 plot(SNR,(mean(ATP3')),'k:>','MarkerSize',10,'LineWidth',3)
 hold off;
 xlabel('SNR (in dB)');
 ylabel('Accuracy (in %)');
 legend('ADTFD','PSNEO','COB');
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(SNR,(mean(AFP1')),'k-o','MarkerSize',10,'LineWidth',3)
 set(gca, 'XDir','reverse','FontSize',25);
 hold on;
 plot(SNR,(mean(AFP2')),'k--s','MarkerSize',10,'LineWidth',3)
 hold on;
 plot(SNR,(mean(AFP3')),'k:>','MarkerSize',10,'LineWidth',3)
 hold off;
 xlabel('SNR (in dB)');
 ylabel('Error (in %)');
 legend('ADTFD','PSNEO','COB');
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 FPR=mean(FPR1'/100);
 TPR=mean(Hitrate1'/100);
 k = length(FPR);
 A = zeros(k,1);

if length(TPR) ~= k
    disp('Length of vectors FPR and TPR must be equal.');
    return;
else
    figure;
    plot(FPR,TPR,'k-','LineWidth',3);
     set(gca,'FontSize',25);
    line([0 FPR(k)], [0 TPR(k)],'color','k');
    line([1 FPR(1)], [1 TPR(1)],'color','k');
    axis([0 1 0 1]);
    
   
    title('ROC Curve');
 xlabel('False Positive Rate(in dB)');
 ylabel('Hit Rate');
end
 for i = 2 : k
    A(i) = (FPR(i - 1) - FPR(i)).*(1 - TPR(i));
end
area1 = 1 - sum(A)

%%%%%%%%%%
 FPR=mean(FPR2'/100);
 TPR=mean(Hitrate2'/100);
 k = length(FPR);
 A = zeros(k,1);

if length(TPR) ~= k
    disp('Length of vectors FPR and TPR must be equal.');
    return;
else
  figure;
    plot(FPR,TPR,'k--','LineWidth',3);
    line([0 FPR(k)], [0 TPR(k)],'color','k');
    line([1 FPR(1)], [1 TPR(1)],'color','k');
    axis([0 1 0 1]);
    title('ROC Curve');
 xlabel('False Positive Rate(in dB)');
 ylabel('Hit Rate');
end

for i = 2 : k
    A(i) = (FPR(i - 1) - FPR(i)).*(1 - TPR(i));
end
area2 = 1 - sum(A) 

%%%%%%%%%%%%%%
 FPR=mean(FPR3'/100);
 TPR=mean(Hitrate3'/100);
 k = length(FPR);
 A = zeros(k,1);

if length(TPR) ~= k
    disp('Length of vectors FPR and TPR must be equal.');
    return;
else
    figure;
    plot(FPR,TPR,'k:','LineWidth',3);
    line([0 FPR(k)], [0 TPR(k)],'color','k');
    line([1 FPR(1)], [1 TPR(1)],'color','k');
    axis([0 1 0 1]);
end
for i = 2 : k
    A(i) = (FPR(i - 1) - FPR(i)).*(1 - TPR(i));
end
area3 = 1 - sum(A) 
 title('ROC Curve');
 xlabel('False Positive Rate(in dB)');
 ylabel('Hit Rate');
 mean(FPR1')
 mean(FPR2')
 mean(FPR3')
