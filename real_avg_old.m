close all;
clear all;
pathD='G:\d\pat006Iktal\';
filenameRoot={'010608ba_0011','010608ba_0023','010608ba_0030'};
maskStart=[	737530, 	757466, 	251198];
maskEnd=[     765163,      768589,      263806];
for f=3:length(filenameRoot)
             for i=3:3
            Ictal=importdata([pathD char(filenameRoot(f)) '_' num2str(i) '.asc']);
%              end
% end
Ictal=interp(decimate(Ictal,2),1);
signal=Ictal(255198/2:257198/2)';
signal=signal-mean(signal);
                signal=signal/norm(signal);
signal=signal(1:512);
% % signal=signal.';
%   x1(i,:)=signal;

             end
end
x=(signal)';
[b,a] = butter(5,0.5/128,'high'); 
x = filtfilt( b,a,x);
% x=diff(abs(x));
% x=[x; 0];
% [b,a] = butter(1,60/128); 
% x = filtfilt( b,a,x);
[tfd,orient]=HTFD_new2_spike(x,3,5,64);
 orient=(orient-1)*45;
  %orient=orient*30;
tfd=imresize(tfd,[length(x) length(x)]);
figure; SetFigDef(16,9,'Times',20);;tfsapl(x,tfd(:,1:4:512),'grayscale','on','sampleFreq','128', 'TFfontSize' , 30);
set(gca, 'FontSize',25);
adtfd1=tfd;
 
 Mask=zeros(size(orient));
 
 Mask( and(orient<100,orient>80))=1;
  %Mask( orient==0)=1;

 Mask(:,1:25)=0;
 Mask(:,end-25:end)=0;
 
 adtfd=Mask.*adtfd1;
 figure; SetFigDef(16,9,'Times',20);; tfsapl(x,adtfd(:,1:4:512),'grayscale','on','sampleFreq','128','TFfontSize' , 30);
set(gca, 'FontSize',25);
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
%   title('(a)');
%  figure;plot(real(x_s),'k-','LineWidth',3);set(gca, 'FontSize',25);
%  xlabel('Sample');
%  ylabel('Amplitude');
%%%%%%%%%%%%%%ADTFD%%%%%%%%%%%%
 figure;plot(real(y),':k','LineWidth',3);set(gca, 'FontSize',25);
 xlabel('Sample');
 ylabel('Amplitude');
 title('(c)');
[p1, kest1]=findpeaks(y,'MINPEAKHEIGHT',0.2) ; 
y= zeros(size(kest1, 1),1); y(kest1)= p1(:);
 figure;plot(real(y),':k','LineWidth',3);set(gca, 'FontSize',25);
 xlabel('Sample');
 ylabel('Amplitude');
%  title('Output of ADTFD after thresholding ');
 %%%%%%%%%%%%%%SNEO%%%%%%%
yy=x(2:end-1).^2-x(3:end).*x(1:end-2);
% yy=filter(gausswin(20,1),1,yy);
yy=yy/max(abs(yy));
% yy=filter(ones(1,5),1,yy);
 yy=[0; yy];
   figure;plot(yy,'k:','LineWidth',3);set(gca, 'FontSize',25);
   xlabel('Sample');
   ylabel('Amplitude');
 title('Output of SNEO before thresholding ');
[p2, kest2]=findpeaks(y,'MINPEAKHEIGHT',0.1) ; 
yy= zeros(size(kest2, 1),1); yy(kest2)= p2(:);
 figure;plot(real(yy),':k','LineWidth',3);set(gca, 'FontSize',25);
 xlabel('Sample');
 ylabel('Amplitude');
 title('Output of SNEO after thresholding ');
 %%%%%%%%%%%%%%%cob%%%%%%%%%%%
   [ft peaktrain]=spbycob1(x, 64);
   peaktrain=peaktrain/max(abs(peaktrain));
   ft=ft/max(abs(ft));
   figure; plot(ft,'k:','LineWidth',3);set(gca, 'FontSize',25);
   xlabel('Sample');
   ylabel('Amplitude'); 
   title('Output of COB before thresholding ');
   figure; plot(peaktrain,'k:','LineWidth',3);set(gca, 'FontSize',25);
   xlabel('Sample');
   ylabel('Amplitude');
 title('Output of COB after thresholding ');