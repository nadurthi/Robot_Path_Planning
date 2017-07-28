function MeasPairs=SensorTask_dummy(MeasPairs,Satellites,Radars,Constants,Tk,Tk1,method)

% MeasPairs{1}(i,j)
%               j=1,rad1  j=2,rad2   j=3,rad3  j=4,rad4   ...
% i=1 sat1
% i=2 sat2
% i=3 sat3
% .
% .
% .


for k=Tk:Tk1
    
%     for i=1:1:Constants.Nsat
%         for j=1:1:Constants.Nrad
%             
%             Satellites{i}
%             Radars{j}
%             
%         end
%     end
    
    M=ones(Constants.Nsat, Constants.Nrad);
    
    MeasPairs{k}=M;




end



