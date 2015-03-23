function scripts = savescripts
%% This script puts the scripts you have in the current folder into variable
% Often, it's useful to save the actual scripts you're running
% with your data. This means you'll be able to go back and see
% what the scripts looked like at the point in time when you
% ran the participant. Only makes files minimally larger.
% Make sure to actually save the output variable with data!
% output is a structure with fields "name" and "content"

all_files = ls('*.m');
num_files = size(all_files, 1);
for i = 1:num_files
    scripts(i).name = all_files(i, :);
    scripts(i).content = fileread(all_files(i, :));
end