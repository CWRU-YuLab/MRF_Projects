function rotMtx = func_rotation(OffRes,b1,nPool,gamma,dt)

nx = 2*pi*gamma*real(b1)*dt;
ny = 2*pi*gamma*imag(b1)*dt;
nz = -2*pi*gamma*OffRes*dt;
rotMtx = zeros(nPool*3,nPool*3);
for ii = 1:nPool
    if(nx==0 && ny==0 && nz(ii)==0)
        rotMtxTemp = eye(3);
    else
        phi = sqrt(nx^2+ny^2+nz(ii)^2);
        hp = phi/2; 
        cp = cos(hp); 
        sp = sin(hp)/phi; 
        ar = cp; 
        ai = -nz(ii)*sp;  
        br = ny*sp; 
        bi = -nx*sp; 
        
        arar = ar*ar;
        aiai = ai*ai;
        arai2 = 2*ar*ai;
        brbr = br*br;
        bibi = bi*bi;
        brbi2 = 2*br*bi; 
        arbi2 = 2*ar*bi; 
        aibr2 = 2*ai*br; 
        arbr2 = 2*ar*br;
        aibi2 = 2*ai*bi;
        
        rotMtxTemp = zeros(3,3);
        rotMtxTemp(1) = arar-aiai-brbr+bibi;
        rotMtxTemp(2) = -arai2-brbi2;
        rotMtxTemp(3) = -arbr2+aibi2;
        rotMtxTemp(4) = arai2-brbi2;
        rotMtxTemp(5) = arar-aiai+brbr-bibi;
        rotMtxTemp(6) = -aibr2-arbi2;
        rotMtxTemp(7) = arbr2+aibi2;
        rotMtxTemp(8) = arbi2-aibr2;
        rotMtxTemp(9) = arar+aiai-brbr-bibi;
    end
    rotMtx((1:3)+3*(ii-1),(1:3)+3*(ii-1)) = rotMtxTemp;
end
