
Wlin=1/1.4;
WoverlBlin=6.3;
kFdslin=1.0;
lMClin=0.0;

Ny=101;
NPhi=151;
NFT=51;

SineApprox=1;

field1 = 'W';
field2 = 'lMC';
field3 = 'lB';
field4 = 'NPhi';
field5 = 'Ny';
field6 = 'PlotFlag';
field7 = 'ShouldSolve';
field8 = 'Gamma';
field9 = 'xZeroCase';
field10 = 'OnCluster';
field12 = 'ChiTop';
field13 = 'ChiBottom';
field14 = 'l';
field15 = 'FirstRun';
field16 = 'L0tinv';
field17 = 'L1';
field18 = 'NyForOutput';

P = struct(field1,0,field2,0,field3,0,field4,0,field5,0,field6,0,field7,0,field8,0,...
    field9,0,field10,0,field12,0,field13,0,field14,0,field15,0,field16,0,...
    field17,0,field18,0);

P.PlotFlag=1;
P.ShouldSolve=1;
P.xZeroCase=0;
P.OnCluster=0;
P.NyForOutput=100;

P.Ny=Ny;
P.NFT=NFT;
P.NPhi=NPhi;


N = length(Wlin);
NlB = length(WoverlBlin);

index=0;

for ilMC =1:length(lMClin)
    lMC=lMClin(ilMC);
    P.lMC=lMC;
    for ikFds=1:length(kFdslin)
        kFds=kFdslin(ikFds);
        P.Gamma = getGamma(NFT,kFds,lMC,SineApprox);
        P.l = 1/P.Gamma(IndexFromN(0,NFT));
        for iW=1:length(Wlin)
            W=Wlin(iW);
            P.W=W;
            for ilB=1:length(WoverlBlin)
                WoverlB = WoverlBlin(ilB);
                lB = W/WoverlB;
                P.lB=lB;
                
                MyResultsTemp = zeros(1,6);
                
                index=index+1;
                
                MyResultsTemp(1,1)=W;
                MyResultsTemp(1,2)=WoverlB;
                MyResultsTemp(1,3)=lMC;
                MyResultsTemp(1,4)=kFds;
                MyResultsTemp(1,5)=Ny;
                MyResultsTemp(1,6)=NPhi;

                result = SolveCompressibleCase(P);

                MyResultsTemp = [MyResultsTemp result];
                fileName=strcat('Results_W_',num2str(W),'_WoverRc_',...
                    num2str(WoverlB),'_kFds_',num2str(kFds),'_lMC_',...
                    num2str(lMC),'_SineApprox_',num2str(SineApprox),'.txt');
                
                dlmwrite(fileName,MyResultsTemp)
                
            end
        end
    end
end

return
