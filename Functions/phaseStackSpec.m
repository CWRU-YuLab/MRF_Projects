function [outputStack] = phaseStackSpec(stack,phase)
%phases all spectra in stack by phase, first order phase corrections applied along the first dimension of stack (unless stack is 1xn)
%phase can be input as either:
%[0th order]
%[0th order, 1st order]
%[0th order, 1st order, 1st order by pivot, pivot point]
%all phases in radians, stack in spectral domain
if size(stack,1)==1
    stack = permute(stack,[2,1,3:ndims(stack)]);
end

outputStack = zeros(size(stack));
for k = 1:size(stack,4)
    for j = 1:size(stack,3)
        for i = 1:size(stack,2)
            outputStack(:,i,j,k) = rephase(stack(:,i,j,k),phase);
        end
    end
end
end