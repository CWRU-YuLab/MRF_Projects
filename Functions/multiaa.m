function out = multiaa(in)

d1 = 1;
dind = 1:5;
for d1 = 1:5
    if size(in,d1)~=1
        break
    end
end
dind(d1) = 1;
dind(1) = d1;
in = permute(in,dind);
out = in;
for i = 1:size(in,2)
    for j = 1:size(in,3)
        for l = 1:size(in,5)
            for k = 1:size(in,4)
                temp = angle(in(:,i,j,k,l)'*ones(size(in,1),1));
                out(:,i,j,k,l)= phaseStackSpec(in(:,i,j,k,l),temp);
            end
        end
    end
end
out = permute(out,dind);