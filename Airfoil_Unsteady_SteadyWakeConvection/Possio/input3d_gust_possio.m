%This script will create the input needed by the Possio code on the SCC
% found in project/aeracous/DOR_POS


% Establish the reduced frequency
k = linspace(0,157,1000);

% Set the Mach number
mach = 0.01;

% Set the marker for transverse gust problem 
trans = 1;

% Set the marker to continue case
count = 1;

fileID = fopen('input3d','w');


for i = 1:length(k)

    fprintf(fileID,'1\n');
    fprintf(fileID,'%f %f\n',mach,k(i));
    fprintf(fileID,'1\n');
    
    
end
fclose(fileID);