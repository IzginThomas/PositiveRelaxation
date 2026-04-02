function P = getProd(S, r)
[N,M] = size(S);
P=zeros(N,N);
for l=1:M
    sl = S(:,l);
    while max(abs(sl)) ~= 0
        slnonzero = sl(sl~=0);
        sbar = min(abs(slnonzero));
        I=find(abs(sl)==sbar);
        i=I(1);
        if sl(i)>0
            [~,J] = min(sl);
            j=J(1);
            P(i,j)=P(i,j) + sl(i)*r(l);
            sl(j) = sl(j) + sl(i);
            sl(i) = 0;
        else
            [~,J] = max(sl);
            j=J(1);
            P(j,i)=P(j,i) + abs(sl(i))*r(l);
            sl(j) = sl(j) + sl(i);
            sl(i) = 0;
        end
    end
end
end