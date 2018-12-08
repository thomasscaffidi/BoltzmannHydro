function result = FixedAcos(x)

% result=acos(x);

if(abs(x-1)<10*eps)
    result=0;
elseif(abs(x+1)<10*eps)
    result=pi;
else
    result=acos(x);
end

end
    