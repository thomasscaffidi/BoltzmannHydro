function result = GetChebyCoord(xmin,xmax,N)

k=N-2:-1:1;
result= [xmin (xmax+xmin)/2+(xmax-xmin)/2*cos((2*k-1)./(2*(N-2))*pi) xmax];




end