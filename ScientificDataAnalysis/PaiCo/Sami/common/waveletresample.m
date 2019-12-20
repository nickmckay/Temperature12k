function y = waveletresample(y,depth, mother)

[C,L] = wavedec(y,depth, mother);
for i = (length(L)-1):-1:2
        C((L(i)+1):L(i+1)) = fastrandperm(C((L(i)+1):L(i+1)));
end

y = waverec(C,L, mother);
