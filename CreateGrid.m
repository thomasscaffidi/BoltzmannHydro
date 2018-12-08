function result = CreateGrid(P)
Ny=P.Ny;
W=P.W;
NPhi=P.NPhi;
lB=P.lB;
OnCluster=P.OnCluster;

%% Create Discretization Grid
YLin=GetChebyCoord(-W/2,W/2,Ny);

k=1:NPhi-2;
PhiZeroLin=[0 pi/2-pi/2*cos((2*k-1)./(2*(NPhi-2))*pi) pi];

phiLim=FixedAcos(1-W/lB);

    if(W<=2*lB)
        
        for i=1:NPhi
            PhiZero=PhiZeroLin(i);
            if(PhiZero<pi-phiLim)
                NewPhiZero = FixedAcos(cos(pi-PhiZero)+W/lB);
                PhiZeroLin = [PhiZeroLin NewPhiZero];
            end
        end
        
        
        if(OnCluster==0)
            PhiZeroLin=uniquetol(PhiZeroLin,10*eps);
        else
             tol=10*eps; %// tolerance
            PhiZeroLin = PhiZeroLin(~any(triu(abs(bsxfun(@minus,PhiZeroLin',PhiZeroLin))<tol,1)));
        end
        NPhi = length(PhiZeroLin);
    
        
        for j=1:length(YLin)
            YTemp=YLin(j);
            phiZeroTemp = FixedAcos( (YTemp+W/2)/lB-1);
            PhiZeroLin = [PhiZeroLin phiZeroTemp];
        end
        if(OnCluster==0)
            PhiZeroLin=uniquetol(PhiZeroLin,10*eps);
        else
            tol=10*eps; %// tolerance
            PhiZeroLin = PhiZeroLin(~any(triu(abs(bsxfun(@minus,PhiZeroLin',PhiZeroLin))<tol,1)));
        end
        %     PhiZeroLin=sort(PhiZeroLin);
        NPhi = length(PhiZeroLin);
        
        %     elseif(W<=4*lB)
    else
        NewN = 2.*round((Ny/2+1)/2)+1;
        
        
        Y1=GetChebyCoord(-W/2,-W/2+2*lB,NewN);
        Y2=GetChebyCoord(-W/2+2*lB,W/2,NewN);
        
        %         Y1=linspace(-W/2,-W/2+2*lB,Ny);
        %         Y2=linspace(-W/2+2*lB,W/2,Ny);
        %         Y3=Y2-2*lB;
        
        YLin = [Y1 Y2 -Y1 -Y2];
        if(OnCluster==0)
            YLin=uniquetol(YLin,10*eps);
        else
            tol=10*eps; %// tolerance
            YLin = YLin(~any(triu(abs(bsxfun(@minus,YLin',YLin))<tol,1)));
        end
        
        %         YLin=uniquetol(YLin,10*eps);
        Ny=length(YLin);
        
        for j=1:length(YLin)
            YTemp=YLin(j);
            if(YTemp>-W/2+2*lB)
                YLin=[YLin YTemp-2*lB -YTemp+2*lB];
            end
            
        end
        
        
        if(OnCluster==0)
            YLin=uniquetol(YLin,10*eps);
        else
            tol=10*eps; %// tolerance
            YLin = YLin(~any(triu(abs(bsxfun(@minus,YLin',YLin))<tol,1)));
        end
        
        Ny=length(YLin);
        
        for j=1:length(YLin)
            YTemp=YLin(j);
            MySol = (YTemp+W/2)/lB-1;
            if(MySol>=-1 && MySol<=1)
                phiZeroTemp = FixedAcos( MySol);
                PhiZeroLin = [PhiZeroLin phiZeroTemp];
            end
            
        end
        %
        if(OnCluster==0)
            PhiZeroLin=uniquetol(PhiZeroLin,10*eps);
        else
            tol=10*eps; %// tolerance
            PhiZeroLin = PhiZeroLin(~any(triu(abs(bsxfun(@minus,PhiZeroLin',PhiZeroLin))<tol,1)));
        end
        NPhi = length(PhiZeroLin);
    end



Points=zeros(1e6,5);
Index=1;
NumberYPoints = zeros(Ny,1);


for iy=1:Ny
    y=YLin(iy);
    PhiZeroStart=1;
    for iPhiZero=PhiZeroStart:NPhi
        phiZero=PhiZeroLin(iPhiZero);
        
        BottomSol = cos(phiZero)-1/lB*(y+W/2);
        
        if(BottomSol<=1+10*eps && BottomSol >= -1-10*eps)
            
            phi = real(FixedAcos(BottomSol));
            if(~(abs(phi-pi)<eps && iy==Ny))
                s=lB*(phi - phiZero);
                
                Points(Index,:)= [ phiZero s phi iy 1];
                Index=Index+1;
                
                NumberYPoints(iy)=NumberYPoints(iy)+1;
                yMax=W/2+lB*(-1-cos(phi));
                
                if(phiZero>pi-phiLim+10*eps    && -phi+2*pi ~= phi && phi~=0) %y<yMax-1000*eps %(y-yMax)/abs(yMax)<-100*eps
                    phi=-phi+2*pi;
                    s=lB*(phi - phiZero);
                    Points(Index,:)= [ phiZero s phi  iy 1];
                    Index=Index+1;
                    
                    
                    NumberYPoints(iy)=NumberYPoints(iy)+1;
                end
            end
            
            
        end
        
        phiZero = pi+phiZero;
        TopSol = cos(phiZero)-1/lB*(y-W/2);
        
        
        if(TopSol<=1+10*eps && TopSol >= -1-10*eps)
            phi = - real(FixedAcos(TopSol))+ 2*pi;
            if(~(abs(phi-2*pi)<10*eps && iy==1))
                s=lB*(phi - phiZero);
                
                Points(Index,:)= [ phiZero s phi iy 2];
                Index=Index+1;
                
                
                NumberYPoints(iy)=NumberYPoints(iy)+1;
                
                
                yMax=-W/2+lB*(1-cos(phi));
                if( phiZero>2*pi-phiLim+10*eps   && phi~=2*pi && phi~=pi) %(y-yMax)/abs(yMax)>100*eps
                    phi = -phi+2*pi;
                    s=lB*(phi - phiZero + 2*pi);
                    
                    Points(Index,:)= [phiZero s phi iy 2];
                    Index=Index+1;
                    
                    
                    NumberYPoints(iy)=NumberYPoints(iy)+1;
                end
            end
            
            
        end
    end
    if(lB<=0.5*W)
        for StartiY=2:Ny
            
            StartY=YLin(StartiY);
            if(StartY<W/2-2*lB)
                Sol = 1-(y-StartY)/lB;
                
                if(Sol<=1+10*eps && Sol >=-1-10*eps)
                    phi=real(FixedAcos(Sol));
                    if(y>-W/2+lB*(1-cos(phi))+eps && y<W/2-lB*(1+cos(phi))-eps)
                        Zone=3;
                        PhiZero=0;
                        s=lB*phi;

                        Points(Index,:)= [StartY s phi iy Zone];
                        Index=Index+1;

                        NumberYPoints(iy)=NumberYPoints(iy)+1;
                        if(phi~=0 && phi~=pi && phi~=2*pi)
                            phi=-phi+2*pi;
                            s=lB*phi;
                            
                            Points(Index,:)= [StartY s phi iy Zone];
                            Index=Index+1;

                            NumberYPoints(iy)=NumberYPoints(iy)+1;
                        end
                    end
                end
            end
        end
    end
end

Size=Index-1;
Points=Points(1:Size,:);

Size = length(Points)

PointsSortedByPhiZeroAndS = sortrows(Points,[5 1 2]);
PointsSortedByPhiZeroAndS = [PointsSortedByPhiZeroAndS (1:Size)'];
PointsSortedByYandPhi = sortrows(PointsSortedByPhiZeroAndS,[4 3]);

Grid=struct('PointsSortedByPhiZeroAndS',PointsSortedByPhiZeroAndS,...
'PointsSortedByYandPhi', PointsSortedByYandPhi,...
'Size',Size,...
'NumberYPoints',NumberYPoints,...
'YLin',YLin,...
'Ny',Ny,...
'NPhi',NPhi);

result = Grid;

if(P.PlotFlag==1)
    Ny
    NPhi
    figure
    scatter(PointsSortedByPhiZeroAndS(:,1)/2/pi,PointsSortedByPhiZeroAndS(:,2))
    grid on
    figure
    colormap jet
    scatter(PointsSortedByYandPhi(:,3)/2/pi,YLin(PointsSortedByYandPhi(:,4)),[],PointsSortedByYandPhi(:,5))
end

% PointsToCheck=Points(:,3:4);
% ToCheckCheck=0;
% if(OnCluster==0)
%     ToCheckCheck=uniquetol(PointsToCheck,eps,'ByRows',true);
% else
%     ToCheckCheck=builtin('_mergesimpts',PointsToCheck,eps*ones(1,2));
% end
% 
% if(length(ToCheckCheck)~=length(PointsToCheck))
%      PointsSortedByYandPhi(:,3:4)
%      length(ToCheckCheck)
%      length(PointsToCheck)
%     error('Twice the same point appears in Points');
% 
% end


end