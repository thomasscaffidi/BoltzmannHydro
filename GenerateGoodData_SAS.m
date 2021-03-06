tic

Npoints=100;
Wlin=logspace(2,-2,Npoints);
% Wlin=0.01;
% Wlin=1;
WoverlBlin=3.9;
kFdslin=0.0;
lMClin=logspace(-2,2,Npoints);
% lMClin=0.05;
UseGuesses=1;
IndexToStart=0;
if isfile('IndexToStart.txt')
     IndexToStart=dlmread('IndexToStart.txt');
else
    dlmwrite('IndexToStart.txt',0);
end

FixeslMCoverW=1;

Ny=75;
NPhi=75;
NFT=1;

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
field19 = 'Guesses';
field20 = 'UseGuesses';
field21 = 'FirstInBatch';

P = struct(field1,0,field2,0,field3,0,field4,0,field5,0,field6,0,field7,0,field8,0,...
    field9,0,field10,0,field12,0,field13,0,field14,0,field15,0,field16,0,...
    field17,0,field18,0,field19,0,field20,0,field21,0);

P.UseGuesses=UseGuesses;
P.PlotFlag=0;
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
AllResults=zeros(N*NlB*length(lMClin)*length(kFdslin),1000);


for ilMC =1:length(lMClin)
    P.FirstInBatch=1;
    for iW=1:length(Wlin)
        
        
        W=Wlin(iW);
        P.W=W;
        lMC=lMClin(ilMC);
        if(FixeslMCoverW==1)
            lMC=lMClin(ilMC)*W;
        end
        P.lMC=lMC;
        for ikFds=1:length(kFdslin)
            kFds=kFdslin(ikFds);
            P.Gamma = getGamma(NFT,kFds,lMC,SineApprox);
            P.l = 1/P.Gamma(IndexFromN(0,NFT));
            
            for ilB=1:length(WoverlBlin)
                WoverlB = WoverlBlin(ilB);
                lB = W/WoverlB;
                P.lB=lB;
                
                MyResultsTemp = zeros(1,6);
                
                index=index+1;
                index
                if(index>=IndexToStart)
                if(index==IndexToStart)
                    P.FirstInBatch=1;
                end
                MyResultsTemp(1,1)=W;
                MyResultsTemp(1,2)=WoverlB;
                MyResultsTemp(1,3)=lMC;
                MyResultsTemp(1,4)=kFds;
                MyResultsTemp(1,5)=Ny;
                MyResultsTemp(1,6)=NPhi;
                
                result = SolveCompressibleCase(P);
                if(P.UseGuesses==1)
                    P.Guesses = result.Guesses;
                end
                
                MyResultsTemp = [MyResultsTemp result.resultArray];
                fileName=strcat('Results_W_',num2str(W),'_WoverRc_',...
                    num2str(WoverlB),'_kFds_',num2str(kFds),'_lMC_',...
                    num2str(lMC),'_SineApprox_',num2str(SineApprox),'.txt');
                
                AllResults(index,1:length(MyResultsTemp))=MyResultsTemp;
                dlmwrite(fileName,MyResultsTemp)
                dlmwrite('IndexToStart.txt',index+1);
                end
                
            end
        end
    P.FirstInBatch=0;    
    end
end
toc
save('AllResults.mat','AllResults')
return
