function dicParam = func_parameter_generation(Variable,nSubDic,nDic)

%check the number of variables
nPool = length(Variable.pool);
nParam = nPool*4+nPool^2; %T1, T2, M0, off-res

%find the values for each parameter, save in cell matrix paramRange
paramRange = cell(1,nParam); %to store the values of all parameters
nDicEntry = 1;
rOffResFlag = false(nPool,nPool); %check if relative off res is used
dMTFlag = false(nPool,nPool); %check if relative MT rate is used
for ii = 1:nPool
    %T1 values of compartment ii
    paramRange{ii} = Variable.pool(ii).T1;
    nDicEntry = nDicEntry*length(paramRange{ii});
    
    %T2 values of compartment ii
    paramRange{ii+nPool} = Variable.pool(ii).T2;
    nDicEntry = nDicEntry*length(paramRange{ii+nPool});
    
    %M0 values of compartment ii
    paramRange{ii+nPool*2} = Variable.pool(ii).M0;
    nDicEntry = nDicEntry*length(paramRange{ii+nPool*2});
    
    %off-resonance values of compartment ii
    if ~isempty(Variable.pool(ii).OffRes) %relative off-res not used
        paramRange{ii+nPool*3} = Variable.pool(ii).OffRes;
        nDicEntry = nDicEntry*length(paramRange{ii+nPool*3});
    else %relative off-res used, #dictionary entry is not changed
        paramRange{ii+nPool*3} = Variable.pool(ii).rOffRes; %assign initial value to be relative freq off-set
        rOffResFlag(ii,Variable.pool(ii).rOffResPool) = true;
    end
    
    %MT rate between compartment ii and compartment jj
    for jj = 1:nPool
        cntMT = sub2ind([nPool nPool],ii,jj); %1D index of (ii,jj) in a matrix with size nPoolxnPool
        if ~isempty(Variable.pool(ii).OutFluxRate{jj})
            if isa(Variable.pool(ii).OutFluxRate{jj},'double')
                paramRange{nPool*4+cntMT} = Variable.pool(ii).OutFluxRate{jj};
                nDicEntry = nDicEntry*length(paramRange{nPool*4+cntMT});
            elseif isa(Variable.pool(ii).OutFluxRate{jj},'char')
                paramRange{nPool*4+cntMT} = 0; %set an intial value
                dMTFlag(ii,jj) = true; %notify that MT(ii,jj) has a rate dependent on MT(jj,ii)*M0(jj)/M0(ii)
            end
        else
            paramRange{nPool*4+cntMT} = 0; %no mt effect
        end
    end
end

%generate all possible combinations of the variables
paramTemp = zeros(nDicEntry,nParam);
%[param1,param2,...paramN] = ndgrid(paramRange{1},paramRange{2},...,paramRange{N})
inputTxt = ' = ndgrid(';
outputTxt = '[';
for ii = 1:(nParam-1)
    inputTxt = strcat(inputTxt,['paramRange{' num2str(ii) '},']);
    outputTxt = strcat(outputTxt,['param' num2str(ii) ',']);
end
inputTxt = strcat(inputTxt,['paramRange{' num2str(nParam) '});']);
outputTxt = strcat(outputTxt,['param' num2str(nParam) ']']);
eval(strcat(outputTxt,inputTxt));
%reshape param1,param2,...,paramN into 1D vector, and combine them in seqDMpar
%format seqDMpar = [param1 param2 ... paramN]
for ii = 1:nParam
    eval(['param' num2str(ii) ' = reshape(param' num2str(ii) ',[],1);']);
    eval(['paramTemp(:,' num2str(ii) ') = param' num2str(ii) ';']);
end

%update the value for relative off-resonance for all entries
[rind,cind] = find(rOffResFlag==true); %compartment r has off-res relative to compartment c
%offRes(r) += offRes(c), offRes(r) was initialized as the freq offset between r and c
paramTemp(:,rind+nPool*3) = paramTemp(:,rind+nPool*3)+paramTemp(:,cind+nPool*3);

%update the relative MT rate for all entries
[rind,cind] = find(dMTFlag==true);
ind1 = sub2ind([nPool,nPool],rind,cind);%1D index of (r,c) in a matrix with size nPoolxnPool
ind2 = sub2ind([nPool,nPool],cind,rind);%1D index of (c,r) in a matrix with size nPoolxnPool
%mt(r,c) = mt(c,r)*m0(c)/m0(r)
paramTemp(:,ind1+nPool*4) = paramTemp(:,ind2+nPool*4).*paramTemp(:,cind+nPool*2)./paramTemp(:,rind+nPool*2);

%select sub dictionary entries
nEntry = ceil(nDicEntry/nSubDic);
if(nDic==nSubDic)
    cEntry = nDicEntry - (nSubDic-1)*nEntry;
else
    cEntry = nEntry;
end
paramTemp = paramTemp((1:cEntry)+nEntry*(nDic-1),:);

dicParam = struct('T1',mat2cell(paramTemp(:,1:nPool),ones(1,cEntry)),...
    'T2',mat2cell(paramTemp(:,(1:nPool)+1*nPool),ones(1,cEntry)),...
    'M0',mat2cell(paramTemp(:,(1:nPool)+2*nPool),ones(1,cEntry)),...
    'OffRes',mat2cell(paramTemp(:,(1:nPool)+3*nPool),ones(1,cEntry)),...
    'MtRate',squeeze(mat2cell(reshape(paramTemp(:,(1:nPool^2)+4*nPool)',nPool,nPool,[]),nPool,nPool,ones(1,cEntry))));