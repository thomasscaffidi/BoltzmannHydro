function MyIntegral = CalculateFT(InputVector,myY,FTindex,Grid,P)

W=P.W;
lB=P.lB;

PointsSortedByYandPhi=Grid.PointsSortedByYandPhi;
NumberYPoints=Grid.NumberYPoints;
YLin=Grid.YLin;
Ny=Grid.Ny;

MyIntegral=0;

SumNy=1;
for iY=1:Ny
    if(iY ~= myY)
        SumNy = SumNy + NumberYPoints(iY);
    else
        y=YLin(iY);
        TempPoints = PointsSortedByYandPhi(SumNy:SumNy+NumberYPoints(iY)-1,:);
        SumNy = SumNy + NumberYPoints(iY);
        PhiLim1=FixedAcos(1-(y+W/2)/lB);
        PhiLim2=-FixedAcos(-1-(y-W/2)/lB)+2*pi;
        
        
        for j=1:NumberYPoints(iY)
            PhiLimTemp1=PhiLim1;
            PhiLimTemp2=PhiLim2;
            
            jPlus=j+1;
            Phi=TempPoints(j,3);
            PhiPlus=0;
            if(j==NumberYPoints(iY))
                %                 dPhi = TempPoints(1,3)+2*pi -  TempPoints(j,3);
                jPlus=1;
                PhiPlus=TempPoints(1,3)+2*pi;
                PhiLimTemp1=PhiLimTemp1+2*pi;
            else
                %                 dPhi = TempPoints(j+1,3) -  TempPoints(j,3);
                PhiPlus = TempPoints(j+1,3);
            end
            
            dPhi=PhiPlus-Phi;
            dPhi1=0.5*dPhi;
            dPhi2=0.5*dPhi;
            alpha1=dPhi1;
            alpha2=dPhi2;
            
            MyIntegral=MyIntegral+alpha1/2/pi*InputVector(TempPoints(j,6))*exp(-1i*FTindex*Phi);
            MyIntegral=MyIntegral+alpha2/2/pi*InputVector(TempPoints(jPlus,6))*exp(-1i*FTindex*PhiPlus);
            
        end
        break
    end
    
end

end