function [] = rb_plotGroups(Controls,Patients)
% Controls is a dataset containing a subject*data matrix of all controls
% Patients is a dataset containing a subject*data matrix of all patients


meanC = mean(Controls);
x = 1:length(meanC);
CImin = meanC-1.96*std(Controls);
CImax = meanC+1.96*std(Controls);

fill([x fliplr(x)],[CImax fliplr(CImin)],'g');
hold on

rb_plotBatch(Patients,1,1);
%set(findobj('Type','line'),'Color','r')

end
