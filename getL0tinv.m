function result = getL0tinv(P,Grid)
lMR=1;
lB=P.lB;
Size=Grid.Size;
PointsSortedByPhiZeroAndS=Grid.PointsSortedByPhiZeroAndS;

UniquePairs = unique(PointsSortedByPhiZeroAndS(:,[1 5]),'rows');
PairIndex=1;
RunningSum=0;
NumberOfPoints=zeros(length(UniquePairs),1);
NumElem=0;
for i=1:Size
    if(i~=Size)
        PhiZero_i = PointsSortedByPhiZeroAndS(i,1);
        PhiZero_ip = PointsSortedByPhiZeroAndS(i+1,1);
        Zone_i = PointsSortedByPhiZeroAndS(i,5);
        Zone_ip = PointsSortedByPhiZeroAndS(i+1,5);
        RunningSum=RunningSum+1;
        if(PhiZero_i ~= PhiZero_ip || Zone_i ~= Zone_ip)
            NumberOfPoints(PairIndex)=RunningSum;
            NumElem=NumElem+RunningSum^2;
            RunningSum=0;
            PairIndex=PairIndex+1;
        end
    else
        RunningSum=RunningSum+1;
        NumElem=NumElem+RunningSum^2;
        NumberOfPoints(PairIndex)=RunningSum;
    end
end

sum(NumberOfPoints)

IndexElem=1;
index=1;
I=zeros(NumElem,1);
J=zeros(NumElem,1);
V=zeros(NumElem,1);

for PairIndex=1:length(UniquePairs)
    Ns = NumberOfPoints(PairIndex);
    LocalPoints = PointsSortedByPhiZeroAndS(index:index+Ns-1,:);
    LocalZone = LocalPoints(1,5);
    if(Ns>1 )
        
        Mat=zeros(Ns,Ns);
        line = zeros(1,Ns);
        
        for lineIndex=2:Ns
            
            s=LocalPoints(lineIndex,2);
            sm=LocalPoints(lineIndex-1,2);
            delta_s=s-sm;
            
            F=zeros(1,Ns);
            F(lineIndex)=1;
            Fm=zeros(1,Ns);
            Fm(lineIndex-1)=1;
            
            line=exp(-(delta_s/lMR))*line...
                +Fm * (lMR/delta_s) * (lMR+exp(-delta_s/lMR)*(-lMR-delta_s)) ...
                + F * (lMR/delta_s) * (-(lMR-delta_s)+lMR*exp(-delta_s/lMR));

%             line=exp(-(delta_s/l))*line+delta_s*0.5*(Fm*exp(-delta_s)+F);
            Mat(lineIndex,:) = line;
        end
        
        if(LocalZone==3)
            s=2*pi*lB;
            sm=LocalPoints(Ns,2);
            delta_s=s-sm;
            
            F=zeros(1,Ns);
            F(1)=1;
            Fm=zeros(1,Ns);
            Fm(Ns)=1;
            
            FinalIntegral=exp(-(delta_s/lMR))*line...
                +Fm * (lMR/delta_s) * (lMR+exp(-delta_s/lMR)*(-lMR-delta_s)) ...
                + F * (lMR/delta_s) * (-(lMR-delta_s)+lMR*exp(-delta_s/lMR));
            
            Mat2 = exp(-LocalPoints(:,2)/lMR)*FinalIntegral/(1-exp(-2*pi*lB/lMR));
            Mat = Mat + Mat2;
            
        end
            

        for ii=1:Ns
            for jj=1:Ns
                V(IndexElem)=Mat(ii,jj);
                I(IndexElem)=LocalPoints(ii,6 );
                J(IndexElem)=LocalPoints(jj,6 );
                IndexElem=IndexElem+1;
            end
        end
    end
    index=index+Ns;
end
I=I(1:IndexElem-1);
J=J(1:IndexElem-1);
V=V(1:IndexElem-1);
result = sparse(I,J,V,Size,Size);

end




