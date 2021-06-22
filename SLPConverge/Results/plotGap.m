%%% scirpt to load each individual file .mat file
load 09
gap=zeros(size(d1,2)-1,1);
fileID = fopen('data_09.txt','w');
fprintf(fileID,'ratio=0.9\n');

for i=3:size(d1,2)
    gap(i-1) = max(abs(d1(i)-d1(i-1))/d1(i-1),abs(d2(i)-d2(i-1))/d2(i-1));
    disp(gap(i-1))
    fprintf(fileID, '%8.4f,%8.4f,%8.4f,%8.4f,%8.4f\n',gap(i-1),d1(i),d2(i),v(6),v(11));
end 
