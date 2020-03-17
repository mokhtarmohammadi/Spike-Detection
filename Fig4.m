
close all;
clear all;
addpath('E:\tfsa_5-5\windows\win64_bin');
n=2047;
 x=cos(2*pi*0.05*(0:255)+2*pi*0.000001*(0:255).^3)+5*cos(2*pi*0.4*(0:255)-2*pi*0.0000001*(0:255).^3);%%+sin(2*pi*0.25*(0:255));%+2*pi*0.001*(0:127).^2)+sin(2*pi*0.4*(0:127));

 %x=50*cos(2*pi*0.05*(0:255)+0*2*pi*0.0000007*(0:255).^3)+100*cos(2*pi*0.4*(0:255)-2*pi*0.0000007*(0:255).^3);%+14*cos(2*pi*0.2*(0:255));
%  x=2*cos(2*pi*0.05*(0:255)+2*pi*0.000001*(0:255).^2);%+2*cos(2*pi*0.4*(0:255));
  x=20*cos(2*pi*0.05*(0:n))+10*cos(2*pi*0.15*(0:n));%+10*cos(2*pi*0.2*(0:255));

 x_s=zeros(1,2048);
 %x_s(40:80:end)=4;
 x_s(40)=1;
 x_s(41)=0;
 x_s(42)=0.0;
 x_s(39)=0;
 x_s(38)=0.0;
 
 
 x_s(120-2:120+2)=x_s(40-2:40+2);
 x_s(400-2:400+2)=x_s(40-2:40+2);
  x_s(600-2:600+2)=x_s(40-2:40+2);
   x_s(1000-2:1000+2)=x_s(40-2:40+2);
    x_s(1200-2:1200+2)=x_s(40-2:40+2);
     x_s(1800-2:1800+2)=x_s(40-2:40+2);
      x_s(2000-2:2000+2)=x_s(40-2:40+2);
 x=x+5*x_s;
 x=awgn(x,30,'measured');
 
 xx=x.';
%x=awgn(x,30,'measured');
  load dataFig4;
  actualpeaks=[40 120 400 600 1000 1200 1800 2000 ];
  fs=1;
 figure;plot(real(x(1:256)),'k-','LineWidth',1);set(gca, 'FontSize',25);
 xlabel('Sample');
 ylabel('Amplitude');
 %title('(a) ');
  axis([0 256 -45 45]);
%  figure;plot(real(x_s),'k-','LineWidth',3);set(gca, 'FontSize',25);
%  xlabel('Sample');
%  ylabel('Amplitude');
%  title('(b) ');
%   axis([0 2048 0 1.1]);

  figure;plot(real(x_s),'k-','LineWidth',3);set(gca, 'FontSize',25);
 xlabel('Sample');
 ylabel('Amplitude');
 title('(b) ');
  axis([0 2048 0 1.1]);

  for seg=1:4
      x=xx((seg-1)*512+1:seg*512);
%  load('shotwithgr');
% ss=e(:,15);
% x=ss(1:4:end);
%  
% x1=denoise(x, 32);
 [adtfd1,orient]=HTFD_new2_spike(x,3,5,100);
 orient=(orient-1)*45;
  %orient=orient*30;
%figure; tfsapl(x,adtfd1(:,1:256),'grayscale','on','sampleFreq','1','Title','(f)', 'TFfontSize' , 30);

  
  
 
 
 Mask=zeros(size(orient));
 
 Mask( and(orient<110,orient>70))=1;
  %Mask( orient==0)=1;

 Mask(:,1:25)=0;
 Mask(:,end-25:end)=0;
 
 adtfd=Mask.*adtfd1;
 for jj=1:length(adtfd1)
 adtfd(jj,:)=adtfd(jj,:)./max(adtfd(jj,:)+0.001);
 end
 
 %tfsapl(x,adtfd)
     %Is=quadtfd(x,length(x)-1,1,'specx',15,'hamm',256);
%figure; tfsapl(x,adtfd)
II=imresize(adtfd,[length(x)/2 length(x)]);
II=adtfd;
%I2=II;
%II(I2<0.2*max(max(II)))
II=II/max(II(:));
II(II>0.02)=1;
II(II<0.02)=0;

MASK=zeros(length(x),length(x));
%MASK(1:end/2,:)=II;
MASK=II;
%T=1;
% mm(mm<T)=0;
% mm(mm>=T)=1;
se = strel('line',5,90);

MASK=imerode(MASK,se);
se = strel('line',5,90);

MASK=imdilate(MASK,se);

mm=sum(MASK);

mm=sum(MASK.*adtfd);

mm=mm/max(mm);

%MASK=imdilate(MASK,se);

%y=ltv_filter((x),MASK);

y = tffilter(MASK,hilbert(x),1:length(x));
y=y.*mm.';
y=mm; 
% 
% y= zeros(size(kest1, 1),1); y(kest1)= p1(:);

%  figure;plot(real(x),'k-','LineWidth',3);set(gca, 'FontSize',25);
%  xlabel('Sample');
%  ylabel('Amplitude');
%  title('(a) ');
% %  axis([0 300 0 1]);
% %  figure;plot(real(x_s),'k-','LineWidth',3);set(gca, 'FontSize',25);
% %  xlabel('Sample');
% %  ylabel('Amplitude');
% %  title('(b) ');
% %   axis([0 2048 0 1.1]);
%%%%%%%%%%%%%%ADTFD%%%%%%%%%%%%
%  figure;plot(real(y),':k','LineWidth',3);set(gca, 'FontSize',25);
%  xlabel('Sample');
%  ylabel('Amplitude');
%  title('Output of ADTFD before thresholding ');
[p1, kest1]=findpeaks(y,'MINPEAKHEIGHT',0.4,'MINPEAKDISTANCE',20) ; 
y= zeros(size(kest1, 1),1); y(kest1)= p1(:);
yy((seg-1)*512+1:seg*512)=[y zeros(1,512-length(y))];
%  hold on;plot(real(y)+(seg-1)*512,':k','LineWidth',3);set(gca, 'FontSize',25);

%  cc=round(actualpeaks*fs);
%     
%      cc=cc(cc<length(x)-10);
    
%     hold on;plot(cc,1.1,'*r','MarkerSize',25,'LineWidth',2)
%        
%     hold off;
%  xlabel('Sample');
%  ylabel('Amplitude');
%  title('(c) ');
  %axis([0 300 0 1.1]);
  end
 figure;plot(real(yy),':k','LineWidth',3);set(gca, 'FontSize',25);

 cc=round(actualpeaks*fs);
    
     cc=cc(cc<length(xx)-10);
  hold on;plot(cc,2,'*r','MarkerSize',25,'LineWidth',2)
  xlabel('Sample');
 ylabel('Amplitude');
 title('(c) ');

  hold off;
   
  %%%%%%%%%%%%%%SNEO%%%%%%%
  clear yy;
yy=xx(2:end-1).^2-xx(3:end).*xx(1:end-2);
yy=filter(gausswin(20,1),1,yy);
yy=yy/max(abs(yy));
% yy=filter(ones(1,5),1,yy);
 yy=[0; yy];
%    figure;plot(yy,'k:','LineWidth',3);set(gca, 'FontSize',25);
%    xlabel('Sample');
%    ylabel('Amplitude');
%  title('Output of SNEO before thresholding ');
[p2, kest2]=findpeaks(yy,'MINPEAKHEIGHT',0.4,'MINPEAKDISTANCE',20) ; 
yy= zeros(size(kest2, 1),1); yy(kest2)= p2(:);
 figure;plot(real(yy),':k','LineWidth',3);set(gca, 'FontSize',25);
 
 hold on;plot(cc,2,'*r','MarkerSize',25,'LineWidth',2);
 %axis([0 300 0 1.1]);
    hold off;

 xlabel('Sample');
 ylabel('Amplitude');
 title('(d) ');
 %%%%%%%%%%%%%%%cob%%%%%%%%%%%
   [ft peaktrain]=spbycob1(xx, 64);
   peaktrain=peaktrain/max(abs(peaktrain));
   ft=ft/max(abs(ft));
%    figure; plot(ft,'k:','LineWidth',3);set(gca, 'FontSize',25);
%    xlabel('Sample');
%    ylabel('Amplitude'); 
%    title('Output of COB before thresholding ');
   [p2, kest2]=findpeaks(peaktrain,'MINPEAKHEIGHT',0.4,'MINPEAKDISTANCE',20); 
yy= zeros(size(kest2, 1),1); yy(kest2)= p2(:);
   figure; plot(yy,'k:','LineWidth',3);set(gca, 'FontSize',25);
 
    hold on;plot(cc,2,'*r','MarkerSize',25,'LineWidth',2);
   % axis([0 300 0 1.1]);
    hold off;
xlabel('Sample');
   ylabel('Amplitude');
 title('(e) ');
%  save Fig2Data x x_s;