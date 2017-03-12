

for j = length(SJ)
    integral(j) = (gamans(j) + gamans(j+1))*SJ(j)/2;
end







% 
% gamma_distance = 0.5*SJ(2:end) + 0.5*SJ(1:end-1);
% gamma_avg = (gamans(2:end) + gamans(1:end-1));
% gamma_integr = gamma_avg.*gamma_distance;
% 
% 
