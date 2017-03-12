

AA = (matA + matB);
BB = matC + matD;

for aa = 1:length(thetpan)+1;
    matE(:,aa) = AA(:,aa).*(-sin(thetpan')) + BB(:,aa).*cos(thetpan');
end

