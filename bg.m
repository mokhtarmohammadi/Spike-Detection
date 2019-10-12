function [Bx,Bc, Px]=bg(data, nfft)


if isempty(data) error('No data to analyze'); end

[Row, Col]= size(data); if (min(Row, Col)~=1) data=data(:); end

len=length(data);
if (floor(len/nfft)==0) error('short data length'); end

%------------------------------- Data Segmentation in M Realization and K Samples each =======

samp=nfft;                % Samples per realization
realz=fix(length(data)/samp);                  % Number of realization
data=data(1:samp*realz);
data=reshape(data,samp,realz);


%------------------------------- Window setting -------------------------------

wind=hann(samp); wind=wind(:);

%------------------------------- Spectrum -------------------------------

if (rem(nfft,2)~=0)
    mrow =fix(nfft/2)+1;
else
    mrow=nfft/2;
end
ncol=mrow;

Px = zeros(nfft,1);
Bx = zeros(mrow,ncol);

mask = hankel([1:mrow],[mrow:mrow+ncol-1] );   % the hankel mask (faster)

Nsamp = [1:samp]';
for Nrealz = 1:realz
    xseg = data(Nsamp)- mean(data(Nsamp)); % Subtract mean value from each record
    xseg   = xseg.*wind;                   % Passing data through window
    Xf     = fft(xseg, nfft)/samp;
    CXf    = conj(Xf);

    Bx  = Bx + (Xf(1:mrow) * Xf(1:ncol).').* CXf(mask);
    Px = Px + Xf.*CXf;

    Nsamp = Nsamp + samp;
end

Px = Px/(realz);        % averaging
Bx = Bx/(realz);        % Bispectrum averaging

%------------------------------- find out the Bicoherence -------------------------------

nm=(Px(1:mrow)*Px(1:ncol).').*Px(mask);
Bc=Bx.^2./(nm);

return
