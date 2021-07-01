function [ output ] = rephase(spectra,phase)
%applies phase correction to spectra
%phase can be input as either:
%[0th order]
%[0th order, 1st order]
%[0th order, 1st order, 1st order by pivot, pivot point]
%all phases in radians, spectra (for first order corrections) must b in frequency domain
phase0 = phase(1);
phase1 = 0;
phase12 = 0;
npts = length(spectra);
pts = [npts-1:-1:0]';
pivot = 0;

if(length(phase)>=2)
    phase1 = phase(2);
    
    if(length(phase)>=4)
        phase12 = phase(3);
        pivot = phase(4);
    end
end
pts12 = [1:npts]'-pivot;
pcorr = phase0 + pts./npts.*phase1+ pts12./npts.*phase12;
specReal = real(spectra);
specImag = imag(spectra);


output = (specReal.*cos(pcorr)-specImag.*sin(pcorr)) + 1i*(specReal.*sin(pcorr)+specImag.*cos(pcorr));
end

