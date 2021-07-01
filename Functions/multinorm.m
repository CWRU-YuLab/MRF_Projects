function [output factor] = multinorm(input)
% if size(input,1)==1
%     %input = input.';
%     input = permute(input,[[2,1],3:ndims(input)]);
% end
% 
% output = zeros(size(input));

d1 = 1;
dind = 1:5;
for d1 = 1:5
    if size(input,d1)~=1
        break
    end
end
dind(d1) = 1;
dind(1) = d1;
input = permute(input,dind);
output = input;

factor = input;

factor(2:end,:,:,:,:,:,:)=[];
for j = 1:size(input,3)
    for i = 1:size(input,2)
        for l = 1:size(input,5)
            for k = 1:size(input,4)
               output(:,i,j,k,l) = input(:,i,j,k,l)./abs(norm(input(:,i,j,k,l)));
               factor(1,i,j,k,l) = abs(norm(input(:,i,j,k,l)));
            end
        end
    end
end
factor = permute(factor,dind);
output = permute(output,dind);