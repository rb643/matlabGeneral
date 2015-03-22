%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to plot matrices with subject by data structure
% Plots all subjects data seperately
%
%% � R.A.I. Bethlehem 2012
%
% Version History:
%
% 01-08-2012: Basic structure and functionality
% 14-08-2012: Randomized line colors
% 15-08-2012: Legend option
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rb_plotBatch(Data,LegendOn,ColourOn)

% check the inputs and sset to defaults if empty
if nargin < 1 || isempty(Data), fprintf('Please enter a valid matrix\n'), return, end;
if nargin < 2 || isempty(LegendOn), LegendOn = 1; end;
if nargin < 3 || isempty(ColourOn), ColourOn = 0; end;

    % get the number of subjects to loop trough
    subjects = size(Data,1);

    % make sure the figure/plot is maintained
    hold on

    % loop through subjects
    for isub = 1:subjects
        h = plot(Data(isub,:));
        
        if ColourOn == 1
        % create a random colour and change the lines color
        color = random('unif',0,1,[1,3]);
        set(h,'Color',color,'LineWidth',1);
        end
        
        % create a legend label for each loop/line
        legendName(isub,:) = sprintf('Subject %03d', isub);
    end
    
    if LegendOn == 1
    % put all the legend labels in the figure
    legend(legendName);
    end
    
    % done plotting so hold off
    hold off
    
end

%% to change all colours after plotting you can also use:
%% set(findobj('Type','line'),'Color','c')
%% where 'c' is the colourcode
