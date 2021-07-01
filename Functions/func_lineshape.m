function lineshape = func_lineshape(type,nSpin,varagin)

if strcmp(type,'Lorentzian')
    %(lw/k)^2/(n^2+(lw/k)^2)
    if isinf(varagin(1))
        lineshape = fftshift([1,zeros(1,nSpin)].');
        return
    end
    fwhm=nSpin.*varagin(1)/2/(nSpin*3);
    
    if mod(nSpin,2)==1
        span = (nSpin-1)/2;
    else
        span = nSpin/2;
    end
    idx = [-span:span].';
    lineshape = (fwhm.^2./(idx.^2+fwhm.^2));
    if mod(nSpin,2)==0
        lineshape(span+1)=[];
    end
    lineshape = reshape(lineshape,[],1);
elseif strcmp(type,'Uniform')
    lineshape = ones(nSpin,1);
end
end