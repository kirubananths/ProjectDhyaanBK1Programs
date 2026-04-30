function idx = encode_state(S, base)

idx = zeros(size(S,1),1);

for j = 1:size(S,2)
    idx = idx + (S(:,j)-1)*base^(j-1);
end

idx = idx + 1;

end