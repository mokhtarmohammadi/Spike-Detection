%% Fig3,4
currentDirectory = pwd;
[upperPath, ~, ~] = fileparts(currentDirectory);
addpath([upperPath '\Core'])
addpath('G:\win64_bin\win64_bin');
%% Signal definition
fs = 32;
n = 0:1/fs:8-1/fs; N = length(n);
sig1 = 2*cos(2*pi*n*0.05*fs);
sig2 = 2*cos(2*pi*n*0.1*fs);
sigm = 0.0020;
B = fir1(16*2,0.08,'high');
s1 = 10*exp(-(fs^2*(n-0.9375).^2)/sigm);
s1 = filter(B,1,s1);
s2 = 10*exp(-(fs^2*(n-2.1875).^2)/sigm);%.*cos(0.1*pi*n);
s2 = filter(B,1,s2);
s3 = 10*exp(-(fs^2*(n- 3.4375).^2)/sigm);%.*cos(0.1*pi*n);
s3 = filter(B,1,s3);
s4 = 10*exp(-(fs^2*(n-4.6875).^2)/sigm);%.*cos(0.1*pi*n);
s4 = filter(B,1,s4);
s5 = 10*exp(-(fs^2*(n- 5.9375).^2)/sigm);%.*cos(0.1*pi*n);
s5 = filter(B,1,s5);
s = s1 + s2 + s3 + s4 + s5;
Signal =  s;
sig=Signal;
Signal = awgn(Signal,0,'measured'); 
load syntethic;
Signal=Signal-mean(Signal);
                Signal=Signal/norm(Signal);
% Signal=abs(filter([1 -1],1,Signal));
% Signal=abs(filter([1 -1],1,Signal));
%   Signal=(exp(abs(Signal))/100000);
% save syntethic  Signal TFD_True;
x=Signal';
% 
%% True TFD
I1 = wvd(s1, N-1) + wvd(s2, N-1) + wvd(s3, N-1) + wvd(s4, N-1) + wvd(s5, N-1);
I2 = wvd(sig1, N-1);
I3 = wvd(sig2, N-1);
TFD_True = I1;

[tfd,orient]=HTFD_new2_spike(x,3,8,54);
 orient=(orient-1)*45;
  %orient=orient*30;
tfd=imresize(tfd,[length(x) length(x)]);
figure; tfsapl(x,tfd(:,1:256),'grayscale','on','sampleFreq','32','Title','(b)', 'TFfontSize' , 30);
% axis([0 300 0 1]);
adtfd1=tfd;
 
 Mask=zeros(size(orient));
 
 Mask( and(orient<115,orient>80))=1;
  %Mask( orient==0)=1;

 Mask(:,1:15)=0;
 Mask(:,end-15:end)=0;
 
 adtfd=Mask.*adtfd1;
 figure; SetFigDef(16,9,'Times',20);tfsapl(x,adtfd(:,1:256),'grayscale','on','sampleFreq','32','Title','(c)', 'TFfontSize' , 30);
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
  title('(a)');
%  figure;plot(real(x_s),'k-','LineWidth',3);set(gca, 'FontSize',25);
%  xlabel('Sample');
%  ylabel('Amplitude');
%%%%%%%%%%%%%%ADTFD%%%%%%%%%%%%
 figure;plot(real(y),':k','LineWidth',3);set(gca, 'FontSize',25);
 xlabel('Sample');
 ylabel('Amplitude');
 title('(d)');
[p1, kest1]=findpeaks(y,'MINPEAKHEIGHT',0.2) ; 
y= zeros(size(kest1, 1),1); y(kest1)= p1(:);
 figure;plot(real(y),':k','LineWidth',3);set(gca, 'FontSize',25);
 xlabel('Sample');
 ylabel('Amplitude');
 title('(e)');
%  title('Output of ADTFD after thresholding ');
 %%%%%%%%%%%%%%SNEO%%%%%%%
yy=x(2:end-1).^2-x(3:end).*x(1:end-2);
% yy=filter(gausswin(20,1),1,yy);
yy=yy/max(abs(yy));
% yy=filter(ones(1,5),1,yy);
 yy=[0; yy];
%    figure;plot(yy,'k:','LineWidth',3);set(gca, 'FontSize',25);
%    xlabel('Sample');
%    ylabel('Amplitude');
%    title('(b)');
%  title('Output of SNEO before thresholding ');
[p2, kest2]=findpeaks(y,'MINPEAKHEIGHT',0.1) ; 
yy= zeros(size(kest2, 1),1); yy(kest2)= p2(:);
 figure;plot(real(yy),':k','LineWidth',3);set(gca, 'FontSize',25);
 xlabel('Sample');
 ylabel('Amplitude');
 title('(f)');
%  title('Output of SNEO after thresholding ');
 %%%%%%%%%%%%%%%cob%%%%%%%%%%%
   [ft peaktrain]=spbycob1(x, 64);
   peaktrain=peaktrain/max(abs(peaktrain));
   ft=ft/max(abs(ft));
% ft(ft<.2)=0;
    [p2, kest2]=findpeaks(peaktrain,'MINPEAKHEIGHT',0.2,'MINPEAKDISTANCE',15); 
yy= zeros(size(kest2, 1),1); yy(kest2)= p2(:);
 figure;plot(real(yy),':k','LineWidth',3);set(gca, 'FontSize',25);

% figure; plot(ft,'k:','LineWidth',3);set(gca, 'FontSize',25);
%    xlabel('Sample');
%    ylabel('Amplitude');
%    title('(b)');
%    title('Output of COB before thresholding ');
%    figure; plot(peaktrain,'k:','LineWidth',3);set(gca, 'FontSize',25);
   xlabel('Sample');
   ylabel('Amplitude');
   title('(g)');
%  title('Output of COB after thresholding ');