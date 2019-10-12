function [ht1, flag]=ifilt(x,nfft)

flag=0;
[bx,bc]=bg(x,nfft);
if sum(real(bx(:)))<0
    [bx,bc]=bg(-x,nfft); flag=1;
end

lb=log(bx); b=ifft(lb);
hx=([b(1,1:end)]); cbhx=hx(2:end)-real(b(1,1)); cbhx=[hx]; cbhx=cbhx-mean(cbhx); chx=exp(cbhx);


% ----------- Fourier Phase Estimation ---------------------

N=nfft/2;

if rem(N,2)==0
    Nby2=N+1;
else
    Nby2=N;
end
r=1;
Nby4=fix(Nby2/2)+1;
R=2*sum(1:Nby4-1)+Nby4-Nby2+1;

psi=zeros(R,1);
amat=zeros(R,Nby2);
for k=2:Nby4
    for l=2:k
        r=r+1;
        psi(r,1)=angle(bx(l,k));
        amat(r,k)=amat(r,k)+1;
        amat(r,l)=amat(r,l)+1;
        amat(r,l+k-1)=amat(r,l+k-1)-1;
    end
end
for k=Nby4+1:Nby2-1
    for l=2:(Nby2+1-k)
        r=r+1;
        psi(r,1)=angle(bx(l,k));
        amat(r,k)=amat(r,k)+1;
        amat(r,l)=amat(r,l)+1;
        amat(r,l+k-1)=amat(r,l+k-1)-1;
    end
end

amat(1,1)=1; psi(1,1)=angle(bx(1,1));

%-------------------- phase unwrapping -----------------------
PHI=[angle(chx(1:end)) 0];
kwrap = fix( (amat*PHI' - psi(1:end)) /(2*pi));
phi = amat(:,1:Nby2-1) \ (psi(1:end) + 2*pi*kwrap);

%--------------- Magnitude & Phase for IFFT ------------------
mag=abs(chx);
mag=[mag]';  mag = [mag(1:Nby2-1); flipud(mag(1:Nby2-2))];
phz = [phi(1:Nby2-1); -flipud(phi(1:Nby2-2))];
chx =( 1./mag) .* exp(-sqrt(-1)*phz);


ht1=ifft(chx);
%ht1=(ifftshift(ht1));
ht1=real(ifftshift(ht1));
if flag==1  ht1=-ht1;  end

return
