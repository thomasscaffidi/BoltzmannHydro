function result = Charac_v8(P)
tic

W=P.W;
lB=P.lB;
NFT=P.NFT;
Gamma=P.Gamma;


Grid = CreateGrid(P);
PointsSortedByPhiZeroAndS=Grid.PointsSortedByPhiZeroAndS;
PointsSortedByYandPhi=Grid.PointsSortedByYandPhi;
NumberYPoints=Grid.NumberYPoints;
Size=Grid.Size;
YLin=Grid.YLin;
Ny=Grid.Ny;


Gamma0 = Gamma(IndexFromN(0,NFT));
Gamma1 = Gamma(IndexFromN(1,NFT));

l=Gamma0^(-1);
lMR=(Gamma0 - Gamma1)^(-1); %should always be one

%% ChiHomZero

ChiHomZero = zeros(Size,1);
for jj=1:Size
    MyData = PointsSortedByPhiZeroAndS(jj,:);
    s=MyData(2);
    phiZero=MyData(1);
    chi=0;
    if(MyData(5)==1)
        chi = (lMR*cos(phiZero)+P.ChiBottom)*exp(-s/lMR);
    elseif(MyData(5)==2)
        chi = (lMR*cos(phiZero)+P.ChiTop)*exp(-s/lMR);
    end
    ChiHomZero(jj) = chi;
end

ChiBulk = -lMR*cos(PointsSortedByPhiZeroAndS(:,3));



%% Generate L0inv and L1
L1=0;
L0tinv=0;
if(P.FirstRun==1)
    L1 = getL1(P,Grid);
    L0tinv = getL0tinv(P,Grid);
else
    L1=P.L1;
    L0tinv=P.L0tinv;
end

% Density1=CalculateFT(ChiHomZero,1,0,Grid,P)
% 
% 
% tt=L1 * ChiHomZero;
% figure
% plot(YLin(PointsSortedByPhiZeroAndS(:,4)),tt,'*')
% 
% error('b')

%% Solve system
  
ChiHom=0;
if(P.ShouldSolve==1)
    tol=1e-5;
    maxit=1000;
    [ChiHom,fl1,rr1,it1,rv1] = bicgstab(@ToSolve,ChiHomZero,tol,maxit,[],[],ChiHomZero);
    fl1
else
    ChiHom = ChiHomZero;
end

resultArray = PostProcess(P,ChiHom,ChiBulk,ChiHomZero,fl1,Grid);

result = struct('resultArray',resultArray,'L0tinv',0,'L1',0);
if(P.FirstRun==1)
    result.L0tinv=L0tinv;
    result.L1=L1;
end


%% tests


% figure
% %     plot(YLin(PointsSortedByPhiZeroAndS(:,4)),temp,'*')
% plot(sort(temp),'*-')
% disp(Ny)
% toc

function OutputVector = ToSolve(InputVector)
    disp('iter')
    temp=L1*InputVector;
% figure
% %     plot(YLin(PointsSortedByPhiZeroAndS(:,4)),temp,'*')
% plot(sort(temp),'*-')
    temp=L0tinv*temp;
    OutputVector = InputVector - temp;
end




end
