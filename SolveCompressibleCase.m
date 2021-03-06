function result = SolveCompressibleCase(P)


if(P.xZeroCase==1)
    P.ChiTop=0;
    P.ChiBottom=0;
    P.FirstRun=1;
    res = Charac_v8(P);
    result = [0 res.resultArray];
else
    
    delta=0.01;
    P.ChiTop=delta;
    P.ChiBottom=-delta;
    P.FirstRun=1;
    resOne = Charac_v8(P);
    yOne = resOne.resultArray(2);
    Guess1 = resOne.ChiHom;
    
    P.L0tinv=resOne.L0tinv;
    P.L1=resOne.L1;
    P.FirstRun=2;
    
    P.ChiTop=-delta;
    P.ChiBottom=delta;
    resMinusOne = Charac_v8(P);
    yMinusOne = resMinusOne.resultArray(2);
    Guess2 = resMinusOne.ChiHom;
    
    b = (yOne - yMinusOne)/2/delta;
    a = yOne - b*delta;
    
    x = -a/b;
    P.ChiTop=x;
    P.ChiBottom=-x;
    P.FirstRun=3;
    res = Charac_v8(P);
    Guess3 = res.ChiHom;
    
    result = struct('resultArray',[x res.resultArray],'Guesses',[Guess1 Guess2 Guess3]);
%     result = [x res.resultArray];
end

end
