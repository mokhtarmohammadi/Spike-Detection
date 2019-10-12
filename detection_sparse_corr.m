close all;
clear all;
fs=200;
M=3;
n=0:127;

addpath('D:\tfsa_5-5\windows\win64_bin');
%well separated components

s1=exp(2*pi*1i*(0.05*n+0.1*n.^2/(2*128)+0.2*n.^3/(128*128*3)));
s2=exp(2*pi*1i*(0.1*n+0.1*n.^2/(2*128)+0.2*n.^3/(128*128*3)));
s3=exp(2*pi*1i*(0.45*n+0.1*n.^2/(2*128)-0.1*n.^3/(128*128*3)));
s4=exp(2*pi*1i*(0.4*n+0.1*n.^2/(2*128)-0.1*n.^3/(128*128*3)));

s = (s1')+ (s2')+ (s3');
s = (s1')+ (s3');

sorig = s;

A=rand(3,1);
X = A*s.';                             % mixed source
% generate noise
SNR=15;
sigma = 10^(-SNR/20);
w = sigma*(randn(M,128) + 1j*(randn(M,128)))/sqrt(2); % noise

X=X+w;
%X=w;
for kkkk=1:M
    r= randsample(128,64) ;
    X(kkkk,r)=0;
end

%Summation of Auto Wigner

I2=quadtfd(X(1,:),length(s)/2-1,1,'wvd',length(X))+quadtfd(X(2,:),length(s)/2-1,1,'wvd',length(X))+quadtfd(X(3,:),length(s)/2-1,1,'wvd',length(X));%+quadtfd(X(4,:),length(s)/2-1,1,'wvd',length(X));
I=quadtfd(X(1,:),length(s)/2-1,1,'mb',0.05,length(X))+quadtfd(X(2,:),length(s)/2-1,1,'mb',0.05,length(X))+quadtfd(X(3,:),length(s)/2-1,1,'mb',0.05,length(X));%+quadtfd(X(4,:),length(s)/2-1,1,'mb',0.05,length(X));

I=post_processing_directional(I2,2,30,60);I(I<0)=0;
%adaptive_optimal_tfd_m;
% Display Wigner distribution
I2(I2<0)=0;
subplot(3,3,1);
imshow(1-I,[])

title('(a)');


% Display Selected TFD






% IF estimation

I9=I/max(max(I));
IF_image_DOA;
%IF_computenew;
K=length(el);


imshow(1-im_bin1',[] )
title('(b)');

% Image segmentation
im_label=zeros(size(I));

for jjj=1:K
    xx=el{jjj};
    for kkk=1:length(xx);
        im_label(xx(kkk,2),xx(kkk,1))=jjj;
    end
    
end
imshow(im_label)
s=conj(s');
dd=[];
I_show=im_label;
% For each segment perform TF filtering to extract components from each
% sensor
% For each segment do the following steps
xx=[];
for jjj=1:K
    IF=zeros(1,length(s));
    IA=IF;
    % Defining TF filter for the given segment
    
    xx=el{jjj};
   
    
    IF(xx(:,1))=(length(s)-xx(:,2))/(2*length(s));
    IA(xx(:,1))=1;
    l=zeros(size(I));
    
    
    Phase=2*pi*filter(1,[1 -1],IF);
    s_dechirp=IA.*exp(1i*Phase);
    
    
    %im_label2=bwmorph(im_label2,'dilate',3);
    
    % For each sensor do the following steps
    
    for iii=1:M
        
        s=(X(iii,:));
        
        %TF filtering for each sensor
        s1 = s.*(s_dechirp);
        s2=fftshift(fft(s1));
        s3=zeros(1,length(s));
        s3(64-2:64+2)=s2(64-2:64+2);
        s1=ifft(ifftshift(s3)).*conj(s_dechirp);
        
        ss(iii,jjj,:)=s1;
        % nly for visualization    TFD of each extracted components
    end
    
    
end

%For each segment compute the covariance matrix and principal eigen vector
sig_den=zeros(M,length(s));
for lll=1:length(ss)
    
    for iii=1:M
        for jjj=1:K
            sig_den(iii,lll)=sig_den(iii,lll)+ss(iii,jjj,lll);
        end
    end
end
%sig_den=X;
AA=cov(sig_den.');

1-det(AA)/prod(diag(AA))

sig_den=X;
AA=cov(sig_den.');

1-det(AA)/prod(diag(AA))