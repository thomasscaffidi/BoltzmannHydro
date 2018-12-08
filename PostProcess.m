function result = PostProcess(P,ChiHom,ChiBulk,ChiHomZero,fl1,Grid)

W=P.W;
lB=P.lB;
Ny=Grid.Ny;
NPhi=Grid.NPhi;
YLin=Grid.YLin;

% YCurrentWithoutL1 =- 2*Jy'*ChiTotalWithoutL1;

XCurrent=zeros(Ny,1);
XCurrentBulk=zeros(Ny,1);
XCurrentHom=zeros(Ny,1);
XCurrentHomZero=zeros(Ny,1);

Density=zeros(Ny,1);
DensityBulk=zeros(Ny,1);
DensityHom=zeros(Ny,1);
DensityHomZero=zeros(Ny,1);

YCurrent=zeros(Ny,1);
YCurrentBulk=zeros(Ny,1);
YCurrentHom=zeros(Ny,1);
YCurrentHomZero=zeros(Ny,1);

for Y=1:Ny
    XCurrent(Y) = - 2*real(CalculateFT(ChiHom,Y,1,Grid,P)) + 1; % 2 factor is for spin
    XCurrentBulk(Y) = - 2*real(CalculateFT(ChiBulk,Y,1,Grid,P)) ;
    XCurrentHom(Y) = - 2*real(CalculateFT(ChiHom,Y,1,Grid,P)) ;
    XCurrentHomZero(Y) = - 2*real(CalculateFT(ChiHomZero,Y,1,Grid,P)) ;
    
    Density(Y) = real(CalculateFT(ChiHom,Y,0,Grid,P));
    DensityBulk(Y) = real(CalculateFT(ChiBulk,Y,0,Grid,P));
    DensityHom(Y) = real(CalculateFT(ChiHom,Y,0,Grid,P));
    DensityHomZero(Y) = real(CalculateFT(ChiHomZero,Y,0,Grid,P));
    
    YCurrent(Y) =  2*imag(CalculateFT(ChiHom,Y,1,Grid,P)) ; % 2 factor is for spin
    YCurrentBulk(Y) =  2*imag(CalculateFT(ChiBulk,Y,1,Grid,P)) ;
    YCurrentHom(Y) =  2*imag(CalculateFT(ChiHom,Y,1,Grid,P)) ;
    YCurrentHomZero(Y) =  2*imag(CalculateFT(ChiHomZero,Y,1,Grid,P)) ;
end

UniformGrid=linspace(-0.5,0.5,P.NyForOutput);
[C,ia,ic] = unique(YLin);
SmoothDensity=interp1(YLin(ia)./W,-(-1/lB*YLin(ia)+Density(ia)'),UniformGrid);
SmoothXCurrent=interp1(YLin(ia)./W,XCurrent(ia)',UniformGrid);
SmoothYCurrent=interp1(YLin(ia)./W,YCurrent(ia)',UniformGrid);

SmoothDensityDeriv=diff(SmoothDensity)./diff(UniformGrid);


if(P.PlotFlag==1)
    figure
    plot(YLin,YCurrent,'*-')
    hold on
    plot(YLin,YCurrentBulk,'*-')
    plot(YLin,YCurrentHomZero,'*-')
    plot(YLin,YCurrentHom,'*-')
    title('Jy')
    legend('Chi','ChiBulk','ChiHomZero','ChiHom')
    
    figure
    plot(YLin,XCurrent,'*-')
    hold on
    plot(YLin,XCurrentBulk,'*-')
    plot(YLin,XCurrentHomZero,'*-')
    plot(YLin,XCurrentHom,'*-')
    title('Jx')
    legend('Chi','ChiBulk','ChiHomZero','ChiHom')
    
    figure
    plot(YLin,Density,'*-')
    hold on
    plot(YLin,DensityBulk,'*-')
    plot(YLin,DensityHomZero,'*-')
    plot(YLin,DensityHom,'*-')
    title('Density')
    legend('Density','DensityBulk','DensityHomZero','DensityHom','Location','BestOutside')
    
    figure
    plot(UniformGrid(1:end-1),SmoothDensityDeriv,'o-')
    title('dRho/dy')
    
end

Xconductivity = trapz(YLin,XCurrent)/W;
Yconductivity = trapz(YLin,YCurrent)/W;
Ey = (Density(end)-Density(1))/W;

result = [Xconductivity Yconductivity Ey P.NyForOutput NPhi UniformGrid ...
    SmoothXCurrent SmoothYCurrent SmoothDensity  SmoothDensityDeriv fl1];
% result=[Xconductivity Yconductivity Ey Ny NPhi YLin XCurrent' YCurrent' Density' fl1];



end