function result=TestApplyL1(P,Grid,InputVector)

NFT=P.NFT;
PointsSortedByYandPhi=Grid.PointsSortedByYandPhi;
NumberYPoints=Grid.NumberYPoints;
Ny=Grid.Ny;
Size=Grid.Size;

NumElem = 0;
for iY=1:Ny
    NumElem=NumElem+NumberYPoints(iY)^2;
end

IndexY=1;
IndexElem=1;
nFT=-NFT:NFT;

I=zeros(NumElem,1);
J=zeros(NumElem,1);
V=zeros(NumElem,1);
result=zeros(Ny,1);

for iY=1:Ny
    LocalPoints = PointsSortedByYandPhi(IndexY:IndexY+NumberYPoints(iY)-1,:);
    LocalPhis = LocalPoints(:,3);
 
    Intervals=diff([LocalPhis; LocalPhis(1)+2*pi]);
    IntegralMask=1/(2*pi)*(0.5*Intervals+0.5*circshift(Intervals,-1));

    U = exp(-1i*nFT'*LocalPhis');

    Mat = ones(NumberYPoints(iY),NumberYPoints(iY))*diag(IntegralMask);
    Mat*InputVector(LocalPoints(:,6))
%     Mat = real(U' * diag(P.Gamma) * U * diag(IntegralMask));

    
    for ii=1:length(LocalPhis)
        for jj=1:length(LocalPhis)
            V(IndexElem)=Mat(ii,jj);
            I(IndexElem)=LocalPoints(ii,6 );
            J(IndexElem)=LocalPoints(jj,6 );
            IndexElem=IndexElem+1;
        end
    end
    
    IndexY = IndexY + NumberYPoints(iY);

end

% result = sparse(I,J,V,Size,Size);


end

