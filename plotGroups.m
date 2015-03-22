function [] = rb_plotGroups(Controls,Patients)
% Controls is a dataset containing a subject*data matrix of all controls
% Patients is a dataset containing a subject*data matrix of all patients

% compute the variables you need
meanC = mean(Controls);
x = 1:length(meanC);
CImin = meanC-1.96*std(Controls);
CImax = meanC+1.96*std(Controls);

% in order to use fill we want to draw a polygon for CI+ in the x direction 
% then go back along the flipped x direction and draw the CI-
fill([x fliplr(x)],[CImax fliplr(CImin)],'g');
hold on

% then we can subsequently use plotbatch to draw the individual patients
rb_plotBatch(Patients,1,1);
%set(findobj('Type','line'),'Color','r')

% NOTE: this will only be somewhat informative if your groups are small
%       otherwise it might be easier to draw another polygon
%       this can be done by using the FaceAlpha transparency option: 
%       fill([x fliplr(x)],[CImax fliplr(CImin)], c, 'FaceAlpha', 0.4)
%       whereby you have to recalculate x and CI for you patient group

end
