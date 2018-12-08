function result = getGamma(NFT,kFds,lMC,SineApprox)

    phi=(0:2*NFT)*2*pi/(2*NFT+1);
    GammaPhi = exp(-4*kFds*sin(phi/2));
    if(SineApprox==1)
        GammaPhi = exp(-2*kFds*abs(wrapToPi_mine(phi)));
    end
    
    
    Gamma = real(fftshift(fft(GammaPhi)));
    
    Gamma0 = Gamma(IndexFromN(0,NFT));
    Gamma1 = Gamma(IndexFromN(1,NFT));
    
    gamma_tr=Gamma0 - Gamma1;
    Gamma=Gamma/gamma_tr;
    
    if(lMC~=0)
        Gamma(IndexFromN(0,NFT))=Gamma(IndexFromN(0,NFT)) +1/lMC;
        Gamma(IndexFromN(1,NFT))=Gamma(IndexFromN(1,NFT)) +1/lMC;
        Gamma(IndexFromN(-1,NFT))=Gamma(IndexFromN(-1,NFT)) +1/lMC;
    end

result=Gamma;
figure
plot(Gamma,'*-')

end