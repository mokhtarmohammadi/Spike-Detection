
close all;
clear all;
 x=cos(2*pi*0.05*(0:255)+2*pi*0.000001*(0:255).^3)+5*cos(2*pi*0.4*(0:255)-2*pi*0.0000001*(0:255).^3);%%+sin(2*pi*0.25*(0:255));%+2*pi*0.001*(0:127).^2)+sin(2*pi*0.4*(0:127));

 x=50*cos(2*pi*0.05*(0:255)+0*2*pi*0.0000007*(0:255).^3)+100*cos(2*pi*0.4*(0:255)-2*pi*0.0000007*(0:255).^3);%+14*cos(2*pi*0.2*(0:255));
%  x=2*cos(2*pi*0.05*(0:255)+2*pi*0.000001*(0:255).^2);%+2*cos(2*pi*0.4*(0:255));
%  x=4*cos(2*pi*0.01*(0:255)+2*pi*0.0001*(0:255).^2)+4*cos(2*pi*0.4*(0:255));

 x_s=zeros(1,256);
 %x_s(40:80:end)=4;
 x_s(40)=1;
 x_s(41)=0;
 x_s(42)=0.0;
 x_s(39)=0;
 x_s(38)=0.0;
 
 
 x_s(120-2:120+2)=x_s(40-2:40+2);
 x_s(200-2:200+2)=x_s(40-2:40+2);
 
 x=x+70*x_s;
 x=x.';
x=awgn(x,15,'measured');
  load Fig2Data1;
  actualpeaks=[40 120 200];
  fs=1;
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

 figure;plot(real(x),'k-','LineWidth',3);set(gca, 'FontSize',25);
 xlabel('Sample');
 ylabel('Amplitude');
 title('(a) ');
%  axis([0 300 0 1]);
 figure;plot(real(x_s),'k-','LineWidth',3);set(gca, 'FontSize',25);
 xlabel('Sample');
 ylabel('Amplitude');
 title('(b) ');
  axis([0 300 0 1.1]);
%%%%%%%%%%%%%%ADTFD%%%%%%%%%%%%
%  figure;plot(real(y),':k','LineWidth',3);set(gca, 'FontSize',25);
%  xlabel('Sample');
%  ylabel('Amplitude');
%  title('Output of ADTFD before thresholding ');
[p1, kest1]=findpeaks(y,'MINPEAKHEIGHT',0.2,'MINPEAKDISTANCE',20) ; 
y= zeros(size(kest1, 1),1); y(kest1)= p1(:);
 figure;plot(real(y),':k','LineWidth',3);set(gca, 'FontSize',25);

 cc=round(actualpeaks*fs);
    
%     cc=cc(cc<length(x)-10);
    
    hold on;plot(cc,1.1,'*r','MarkerSize',25,'LineWidth',2)
       
    hold off;
 xlabel('Sample');
 ylabel('Amplitude');
 title('(c) ');
  axis([0 300 0 1.1]);
 %%%%%%%%%%%%%%SNEO%%%%%%%
yy=x(2:end-1).^2-x(3:end).*x(1:end-2);
yy=filter(gausswin(20,1),1,yy);
yy=yy/max(abs(yy));
% yy=filter(ones(1,5),1,yy);
 yy=[0; yy];
%    figure;plot(yy,'k:','LineWidth',3);set(gca, 'FontSize',25);
%    xlabel('Sample');
%    ylabel('Amplitude');
%  title('Output of SNEO before thresholding ');
[p2, kest2]=findpeaks(yy,'MINPEAKHEIGHT',0.2,'MINPEAKDISTANCE',20) ; 
yy= zeros(size(kest2, 1),1); yy(kest2)= p2(:);
 figure;plot(real(yy),':k','LineWidth',3);set(gca, 'FontSize',25);
 
 hold on;plot(cc,1.1,'*r','MarkerSize',25,'LineWidth',2);
 axis([0 300 0 1.1]);
    hold off;

 xlabel('Sample');
 ylabel('Amplitude');
 title('(d) ');
 %%%%%%%%%%%%%%%cob%%%%%%%%%%%
   [ft peaktrain]=spbycob1(x, 64);
   peaktrain=peaktrain/max(abs(peaktrain));
   ft=ft/max(abs(ft));
%    figure; plot(ft,'k:','LineWidth',3);set(gca, 'FontSize',25);
%    xlabel('Sample');
%    ylabel('Amplitude'); 
%    title('Output of COB before thresholding ');
   [p2, kest2]=findpeaks(peaktrain,'MINPEAKHEIGHT',0.2,'MINPEAKDISTANCE',20); 
yy= zeros(size(kest2, 1),1); yy(kest2)= p2(:);
   figure; plot(yy,'k:','LineWidth',3);set(gca, 'FontSize',25);
 
    hold on;plot(cc,1.1,'*r','MarkerSize',25,'LineWidth',2);
    axis([0 300 0 1.1]);
    hold off;
xlabel('Sample');
   ylabel('Amplitude');
 title('(e) ');
%  save Fig2Data x x_s;