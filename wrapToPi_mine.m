function result = wrapToPi_mine(lambda)
result=lambda;
for i=1:length(lambda)
if(lambda(i)>pi)
    result(i)=lambda(i)-2*pi;
end
end
end
