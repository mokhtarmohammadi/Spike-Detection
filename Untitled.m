clear all;
close all;
n=0:255;
addpath('D:\tfsa_5-5\windows\win64_bin');
pathD='G:\d\pat002Iktal\';
filenameRoot={'010403aa_0015','010403aa_0018','010403aa_0021'};
maskStart=[819373	, 252098, 767172];
maskEnd=[    857094,  269660, 802660];


for f=2:2%length(filenameRoot)
             for i=1:6
            Ictal=importdata([pathD char(filenameRoot(f)) '_' num2str(i) '.asc']);
             

% Ictal = exp((((diff(Ictal))))/100000);
Ictal=resample(Ictal,128,256);
x=Ictal(255598/2:255598/2+511)';
s1=cos(2*pi*(0.05*n+0.05*n.^2/(2*128)+0.05*n.^3/(128*128*3)));
s2=cos(2*pi*(0.15*n+0.05*n.^2/(2*128)+0.05*n.^3/(128*128*3)));
s3=cos(2*pi*(0.45*n));
s=s1+s2;%+s3;
r=rand(1,256);
r(r<0.5)=0;
r(r>=0.5)=1;
s=s.*r;
s=x;
[Inew1,Iorient]=HTFD_new2(s,3,8,84);
% figure; tfsapl(s,Inew1)
figure; tfsapl(s,Inew1(:,1:4:512),'grayscale','on','sampleFreq','128','Title','(ch)', 'TFfontSize' , 20);

%[Inew2,Iorient]=HTFD_new2(s,2,30,156);
%[Inew1,Iorient]=HTFD_new2(s,2,30,196);
%Inew1=min(min(Inew1,Inew2),Inew3);
Iorient=Iorient*3;
  x=s.';
  % adaptive_optimal_tfd;
%Inew=I_max_new/sum(sum((I_max_new)));
%Inew1(and(70<Iorient, Iorient<110))=0;
%Inew1 = quadtfd(s, length(s)/4-1, 1, 'mb',0.05,length(s));
I_max_new=Inew1;
IF_COMPUTE_new_edge_link_1;
Inew1=Inew1.*ei';
% figure;tfsapl(s,Inew1)
figure; tfsapl(s,Inew1(:,1:4:512),'grayscale','on','sampleFreq','128','Title','(ch)', 'TFfontSize' , 20);
%figure;
%tfsapl(s,Inew)
             end
end