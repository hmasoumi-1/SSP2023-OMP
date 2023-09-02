function A = csmat(Pvec,UnitDFTmtx,CircShifts)

Pvec = reshape(Pvec,[1,numel(Pvec)]); % row vector
P = [];
for itr = 1:length(CircShifts)
    tmp = circshift(Pvec,CircShifts(itr));
    P = [P;tmp];
end
A = P*UnitDFTmtx;

end