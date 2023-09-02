function mu = mutual_coherence(A)

An = normalize(A,'norm',2);
Gram = (An')*An;
Gram_zero_diag = Gram - diag(diag(Gram));
mu = max(max(abs(Gram_zero_diag)));

end