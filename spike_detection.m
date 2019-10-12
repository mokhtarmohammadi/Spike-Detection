
close all;
clear all;
addpath('D:\classes\time-frequency\tfsa_6-0\tfsa_6-0\windows\win64_bin');
SNR=[-15 -10 -5 0 5 10 15];
nfft=32;
for nn=1:7
    if nn==1
        
        [finalsignal actualpeaks stimes_1 rngreturnstatus] = makenoisysamples();
    else
        [finalsignal actualpeaks stimes_1 rngreturnstatus] = makenoisysamples('ReuseTimes',stimes_1,'ReuseRNGstate',rngreturnstatus,'NoiseSNR',SNR(nn));
    end
    
    % x1=finalsignal;
    finalsignal=finalsignal-mean(finalsignal);   finalsignal=finalsignal./max(abs(finalsignal)); finalsignal =finalsignal(:);
    x=finalsignal(1:1024);
    %  x=cos(2*pi*0.05*(0:255)+2*pi*0.000001*(0:255).^3)+5*cos(2*pi*0.4*(0:255)-2*pi*0.0000001*(0:255).^3);%%+sin(2*pi*0.25*(0:255));%+2*pi*0.001*(0:127).^2)+sin(2*pi*0.4*(0:127));
    %
    %  x=2*cos(2*pi*0.05*(0:255)+2*pi*0.000001*(0:255).^2)+2*cos(2*pi*0.4*(0:255));
    %  x_s=zeros(1,256);
    %  x_s(40:80:end)=4;
    %  x=x+x_s;
    %  x=x.';
    %  x=signals.';
    %  x1=decimate(X,10);
    %  x=x1(1:512);
    %x=awgn(x,15,'measured');
    
    %  [adtfd1,orient]=HTFD_new2_spike(x,2,15,94);
    %  adtfd1=optshrink(adtfd1,2);
    %  orient=orient*3;
    denoised_signal=denoise(finalsignal, nfft);
    if SNR(nn) <25
        [adtfd1,orient]=HTFD_new2_spike(denoised_signal,3,5,100);
         orient=(orient-1)*45;
    else
    [adtfd1,orient]=HTFD_new2_long(denoised_signal,3,20,256);
    orient=(orient)*15;
    end
%     if (SNR(nn)) <10
%         [adtfd1,relmse_hat,mse_hat] = optshrink((adtfd1),32);
%     end
%     
   
    for jj=1:length(adtfd1)
        %adtfd1(jj,:)=adtfd1(jj,:)./max(adtfd1(jj,:)+0.01);
    end
    
    
    Mask=zeros(size(orient));
    
    Mask( and(orient<110,orient>70))=1;
    %Mask( orient==0)=1;
    
    Mask(:,1:15)=0;
    Mask(:,end-15:end)=0;
    
    adtfd=Mask.*adtfd1;
    %tfsapl(x,adtfd)
    %Is=quadtfd(x,length(x)-1,1,'specx',15,'hamm',256);
    %figure; tfsapl(x,adtfd)
    II=imresize(adtfd,[length(x)/2 length(x)]);
    %I2=II;
    %II(I2<0.2*max(max(II)))
    II=II/max(II(:));
    II(II>0.05)=1;
    II(II<0.05)=0;
    
    MASK=zeros(length(x),length(x));
    MASK(1:end/2,:)=II;
    mm=sum(MASK);
    mm(mm<30)=0;
    mm(mm>=30)=1;
    se = strel('line',20,90);
    MASK=bwmorph(MASK,'clean');
    MASK=bwmorph(MASK,'dilate');
    
    %MASK=imdilate(MASK,se);
    
    %y=ltv_filter((x),MASK);
    
    y = tffilter(MASK,hilbert(x),1:length(x));
    y=y.*mm.';
    detected=real(y);
    detected=detected/max(abs(detected));
    %%%%%%%%%%%%%%%%%%%
    if (SNR(nn)) <10
    detected=find(detected>0.001*max(abs(detected)) );
%     detected=find(detected>4*mean(abs(detected)) );
    else
      detected=find(detected>0.001*max(abs(detected)) );  
%       detected=find(detected>0.5*mean(abs(detected)) );
    end
    if ~isempty(detected); c_detected = zeros(length(finalsignal), 1); c_detected (detected)= 1;end;
    kest1=find(c_detected(1:end)>0); 
    figure;plot(kest1,1,'*r');
    title('ADTFD','FontSize',20);
    hold on;
    plot(finalsignal);set(gca,'FontSize',20);
    hold on;plot(actualpeaks*20480,-1,'Og')
    hold off;
    
    %%%%%%%%%%%%%%%%%%%%%%%  SNEO  %%%%%%%%%%%
    yy=x(2:end-1).^2-x(3:end).*x(1:end-2);
    yy=filter(gausswin(20,1),1,yy);
    yy=[0;0; yy]';
    yy=yy/max(abs(yy));
    yy=find(yy>0.6*max(abs(yy)) );
    if ~isempty(yy); c_yy = zeros(length(finalsignal), 1);c_yy(yy)= 1;end;
     kest2=find(c_yy(1:end)>0); 
    figure; plot(finalsignal);set(gca,'FontSize',20);
    hold on;
    plot(kest2,1,'Xr');
    %
    title('SNEO');
    hold on;plot(actualpeaks*20480,-1,'Og');
    hold off;
    %%%%%%%%%%%%%%%%%    COB    %%%%%%%
    peaktrain=spbycob(x, 128);
    %=========   thresholding ===========
    peaktrain=peaktrain/max(abs(peaktrain));
    thrl=0.3*max(peaktrain);
    spikes= find(peaktrain>thrl);
    if ~isempty(spikes); c_spike = zeros(length(finalsignal), 1); c_spike (spikes)= 1; end;
        kest3=find(c_spike(1:end)>0); 
%     id=find(diff(kest1)< (mintime_m*fs));
    %================Estimate Spikes, Missing & Insertion===========
    
    %c_spike =findimpulse(dataset, 256, thrl); kest1=[find(c_spike(124:end)>0)]; id=find(diff(kest1)< (mintime_m*fs));
    %c_spike =neuron_impulse(dataset, 256, thrl);kest1=[find(c_spike(124:end)>0)]; id=find(diff(kest1)< (mintime_m*fs));
    % Estimated Position of Peaks
    figure;plot(finalsignal);set(gca,'FontSize',20);
    title('COB');
    hold on;
    plot(kest3,1,'*r');
    hold on;plot(actualpeaks*20480,-1,'Og')
    hold off;
    clear detected kest1 kest2 kest3;
end