function [denoised_signal]=denoise(signal, nfft)

signal=signal-mean(signal);   signal=signal./max(abs(signal)); signal =signal(:);
lt=length(signal);
% ------------------ Inverse Filter -----------------------------
[iht]=ifilt(signal,nfft);

%----------------------Noisy Impulse find -----------------------
signal=[signal(1:end-1)+signal(2:end); 0]*.5;
ft=conv(signal,iht); 
ft=ft(:); ft=[zeros(nfft/2+6,1); ft((nfft+1):(end-nfft))]; 

%================Denoising =====================
xtra0=length(ft)-nfft*(floor(length(ft)/nfft)); 
if xtra0>0; ft=[ft; zeros(2^nextpow2(xtra0)-xtra0,1)]; end;    % Zero padding to use wevelet function

%[s] = wden(ft,'rigrsure','s','mln',n,'bior1.5');

temp=ft; [a,d]=swt(ft, 3, 'coif1');  
if abs(skewness(temp))<1
    for n=1:3
        s=iswt(a(n,:), zeros(1, length(d)), 'coif1')';
        if abs(skewness(s))> abs(skewness(temp))  
            temp=s;
            if abs(skewness(s))>=1; break; end
        end
    end
end

if abs(skewness(temp))<1
    for n=1:3
        s=iswt( zeros(1, length(d)), d(n,:), 'coif1')';
        if abs(skewness(s))> abs(skewness(temp))    
            temp=s; %disp('used D');
            if abs(skewness(s))>=1; break; end
        end
    end
end
ft=temp; 

%----------------------Noise Suppressing -------------
f(lt,1)=0;  f(1:length(ft))=ft;
t1=f.*(f>=0);                    %t1 matches at s(nfft/2+7:end -nfft/2-6) for WT &   s(nfft/2:end -nfft/2)
t2=t1.*t1;
t3=t2.*t1;g=sort(t3,'descend'); t=(10*t3/min(g(1:10)));
denoised_signal=t3;
return

