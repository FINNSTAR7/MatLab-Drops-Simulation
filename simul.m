%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINNER
% September 15th, 2020
% Drop Simulation

% For a specified number of iterations and desired probibilities simul()
% will count the number of times it takes to succeed for each probability
% at least once, plot the data, and calculated the average/expected value
% as well as the values at which 99%, 99.9%, and 99.99% certainty is
% reached.

% Inputs:
% Sample Size --- Number of players to simulate.
% Probabilities - The probability of each desired part.
%                 (separated by spaces, not commas)
% Parts --------- The number of parts desired per probability.
%                 (separated by spaces, not commas)

% Outputs:
% None

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -Error Checks- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defs = {}; n = 0;
if exist('samples', 'var')
    defs{1} = num2str(samples);
else
    defs{1} = '500000';
end
if exist('probabilities', 'var')
    defs{2} = num2str(probabilities);
    n = length(probabilities);
else
    defs{2} = '';
end
if exist('parts', 'var')
    defs{3} = num2str(parts);
elseif n > 0
    defs{3} = num2str(ones([1 n]));
else
    defs{3} = '';
end

inputs = inputdlg({'Number of Samples',...
    'Probabilities', 'Parts'},...
    'Sample Size', [1 50; 1 50; 1 50], defs);

if isempty(inputs)
    clearvars -except samples probabilities parts
    return
end

if isempty(inputs{1})
    if ~exist('samples', 'var')
        samples = 500000;
    end
else
    samples = str2num(inputs{1});
end
if floor(samples) ~= ceil(samples) || samples <= 0
    clearvars -except samples probabilities parts
    error('Sample size must be a positive whole number.');
end
assignin('base', 'samples', samples);

if isempty(inputs{2})
    if ~exist('probabilities', 'var')
        clearvars -except samples probabilities parts
        error('No probability array was given nor already exists.');
    end
else
    probabilities = str2num(inputs{2});
end
n = length(probabilities);
for i=1:n
    if(probabilities(i) <= 0)
        clearvars -except samples probabilities parts
        error('Probabilities must be greater than zero.');
    end
end
assignin('base', 'probabilities', probabilities); probs = probabilities;

if isempty(inputs{3})
    if ~exist('parts', 'var') || length(parts) ~= n
        parts = ones([1 n]);
    end
else
    parts = str2num(inputs{3});
end
for i=1:length(parts)
    if floor(parts(i)) ~= ceil(parts(i)) || parts(i) <= 0
        clearvars -except samples probabilities parts
        error('Part counts must be positive whole numbers.');
    end
end
assignin('base', 'parts', parts);
clear defs inputs

%%%%%%%%%%%%%%%%%%%%%%%% -Begin Simulation & Plot- %%%%%%%%%%%%%%%%%%%%%%%%
if sum(probs) >= 1 % if total probability exeeds 100%, use their weights
    probs = probs / sum(probs);
end
null_prob = 1-sum(probs);

if  null_prob < 0
    null_prob = 0;
end

currCount = zeros([samples n+1]); y_vals = [];
num = 0; bool = false;

tic; t = 0;
str = sprintf(['Iteration: 1',...
    '\nSuccesses: ' addComma(num) '/' addComma(samples), '\nProgress: 0%%'...
    '\nTime Elapsed: 00 Seconds, 00 Minutes, 00 Hours, 0 Days']);
f = waitbar(0, str, 'Name', 'Simulating...', 'CreateCancelBtn',...
    'setappdata(gcbf, ''canceling'', 1)');
setappdata(f, 'canceling', 0); k = 1;
while (num < samples) && (k < 1000000)
    if getappdata(f, 'canceling')
        delete(f);
        clearvars -except samples probabilities parts
        return
    end
    
    vals = randsample(0:n, samples-num, true, [null_prob probs])';
    
    for i=1:n+1
        if getappdata(f, 'canceling')
            delete(f);
            clearvars -except samples probabilities parts
            return
        end
        
        currCount(:,i) = currCount(:,i) + (vals == i-1);
        t = updateProgress(num, samples, t, f, k);
    end
    
    rows = []; prevArr = [];
    for i=1:n
        if getappdata(f, 'canceling')
            delete(f);
            clearvars -except samples probabilities parts
            return
        end
        
        [arr, ~] = find(vals == i);
        rows(length(prevArr)+1:length(prevArr)+length(arr),1) = arr;
        prevArr = arr;
        t = updateProgress(num, samples, t, f, k);
    end
    
    deleteRows = zeros([1 length(rows)]);
    for i=1:length(rows)
        if getappdata(f, 'canceling')
            delete(f);
            clearvars -except samples probabilities parts
            return
        end
        
        for j=2:n+1
            if currCount(rows(i),j) >= parts(j-1)
                bool = true;
            else
                bool = false;
                break;
            end
            
            t = updateProgress(num, samples, t, f, k);
        end
        
        if bool
            index = sum(currCount(rows(i),:));
            if index > length(y_vals)
                y_vals(index,1) = 0;
            end
            y_vals(index,1) = y_vals(index,1) + 1;
            deleteRows(i) = rows(i);
            
            num = num + 1;
            bool = false;
        end
        
        t = updateProgress(num, samples, t, f, k);
    end
    
    deleteRows(deleteRows == 0) = [];
    currCount(deleteRows,:) = [];
    k = k + 1;
end
waitbar(1, f, 'Simulation Finished!'); pause(1); delete(f); k = k - 1;
clearvars -except samples probabilities parts y_vals k n probs

x_vals = (sum(parts):k)';
if sum(parts) > 1
    y_vals(1:sum(parts)-1) = [];
end

if n == sum(parts)
    average = averageRun(probs);
else
    average = sum(x_vals.*y_vals)/samples;
end

mode = x_vals(y_vals == max(y_vals));
[average_p, mode_p] = trapezoid_data(average, mode, x_vals, y_vals);

[median, ng] = reverse_trapezoid(x_vals, y_vals);
range = round((ng(3) - ng(1))/2, 0);

varience = sum(y_vals.*(x_vals.^2))/samples - average^2;
std_err = sqrt(varience);

d_x_vals = x_vals(y_vals > 0); d_y_vals = y_vals(y_vals > 0);
x_vals = linspace(d_x_vals(1), d_x_vals(end), 10000); y_vals = x_vals;
[~, i] = max(x_vals >= ng(1));
y_vals1 = interp1(d_x_vals, d_y_vals, x_vals, 'spline');
y_vals2 = interp1(d_x_vals, d_y_vals, x_vals, 'pchip');
y_vals(1:i) = y_vals1(1:i); y_vals(i+1:end) = y_vals2(i+1:end);
y_vals = abs(y_vals);

time = toc/60/60/24; days = fix(time); hrs = fix((time - days)*24);
mins = fix(((time - days)*24 - hrs)*60);
secs = (((time - days)*24 - hrs)*60 - mins)*60;

%
figure; % from here on is mostly just plotting customization
cmap = colormap(jet(length(x_vals)));
set(gcf, 'position', [10 50 800 600]); d_y_vals = d_y_vals/samples;
b = bar(d_x_vals, d_y_vals, 'FaceColor', 'flat'); hold on
c = linspace(0, 1, length(x_vals)); y_vals = y_vals/samples;
scatter(x_vals, y_vals, 0.5, c); c = interp1(x_vals, c, d_x_vals);
b.CData = c';
yt = get(gca, 'YTick'); delta = (yt(end) - yt(1))/(length(yt) - 1);
yLim = [0 yt(end)+2*delta]; set(gca, 'YLim', yLim);

val = mode; [~, i] = max(x_vals >= val);
plot([val val], yLim, 'Color', cmap(i, :));
text(val, yLim(2)*0.925, sprintf([' Mode'...
    '\n x = ' num2str(val)...
    '\n y = ' num2str(max(y_vals)*100,'%.3f')...
    '%%\n A = ' num2str(mode_p*100,'%.3f') '%%']));

val = median; [~, i] = max(x_vals >= val);
plot([val val], yLim, 'Color', cmap(i, :));
text(val, yLim(2)*0.775, sprintf([' Median'...
    '\n x = ' num2str(val)...
    '\n y = ' num2str(freq(val, x_vals, y_vals)*100,'%.3f')...
    '%%\n A = 50%%']));

val = average; [~, i] = max(x_vals >= val);
plot([val val], yLim, 'Color', cmap(i, :));
text(val, yLim(2)*0.625, sprintf([' Average'...
    '\n x = ' num2str(val)...
    '\n y = ' num2str(freq(val, x_vals, y_vals)*100,'%.3f')...
    '%%\n A = ' num2str(average_p*100,'%.3f') '%%']));

val = ng(1); [~, i] = max(x_vals >= val);
plot([val val], yLim, 'Color', cmap(i, :));
text(val, yLim(2)*0.925,...
    sprintf([' x = ' num2str(val)...
    '\n y = ' num2str(freq(val, x_vals, y_vals)*100,'%.3f')...
    '%%\n A = 99%%']));

val = ng(2); [~, i] = max(x_vals >= val);
plot([val val], yLim, 'Color', cmap(i, :));
text(val, yLim(2)*0.775,...
    sprintf([' x = ' num2str(val)...
    '\n y = ' num2str(freq(val, x_vals, y_vals)*100,'%.3f')...
    '%%\n A = 99.9%%']));

val = ng(3); [~, i] = max(x_vals >= val);
plot([val val], yLim, 'Color', cmap(i, :));
text(val, yLim(2)*0.625,...
    sprintf([' x = ' num2str(val)...
    '\n y = ' num2str(freq(val, x_vals, y_vals)*100,'%.3f')...
    '%%\n A = 99.99%%']));

title(sprintf(['Number of Runs for ' addComma(samples) ' Samples']));
xlabel('Number of Runs'); ylabel('Frequency');
xlim([sum(parts)-1 k]); ax = gca; ax.YGrid = 'on'; ax.XGrid = 'on';
%set(gca, 'xscale', 'log');

str = ['  CI_{50%%}= ' char(177), num2str(0.6745*std_err,'%.5f')...
    '    CI_{75%%}= ' char(177), num2str(1.15*std_err,'%.5f')...
    '    CI_{90%%}= ' char(177), num2str(1.645*std_err,'%.5f')...
    '    CI_{95%%}= ' char(177), num2str(1.96*std_err,'%.5f')...
    '    CI_{99%%}= ' char(177), num2str(2.576*std_err,'%.5f')];
text(x_vals(1), yLim(2)*0.0275, sprintf(str));
hold off
%}

fprintf(['   Time Elapsed      = ' num2str(secs) ' Seconds, '...
    num2str(mins) ' Minutes, ' num2str(hrs) ' Hours, ' num2str(days)...
    ' Days\n']);
fprintf(['   Expected          = ' num2str(average) ' - '...
    num2str(average_p*100) '%%\n']);
fprintf(['   Median            = ' num2str(median) ' - 50%%\n']);
fprintf(['   Mode              = ' num2str(mode) ' - '...
    num2str(mode_p*100) '%%\n']);
fprintf(['   Nearly Guaranteed = ' num2str(round(ng(2))) ' '...
    char(177) ' ' num2str(range) '\n\n']);
clearvars -except samples probabilities parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t] = updateProgress(num, samples, t, f, k, del)
if(nargin < 7)
    del = 1;
end
t_curr = toc;
if (t_curr - t) >= del
    m = 100*num/samples;
    time = t_curr/60/60/24; days = fix(time); hrs = fix((time - days)*24);
    mins = fix(((time - days)*24 - hrs)*60);
    secs = fix((((time - days)*24 - hrs)*60 - mins)*60);
    te = sprintf('%02d Seconds,', secs);
    te = strcat(te, sprintf(' %02d Minutes,', mins));
    te = strcat(te, sprintf(' %02d Hours,', hrs));
    te = strcat(te, sprintf(' %0d Days', days));
    str = sprintf(['Iteration: ' addComma(k),...
        '\nSuccesses: ' addComma(num) '/' addComma(samples),...
        '\nProgress: ' num2str(floor(m)) '%%', '\nTime Elapsed: ' te]);
    waitbar(m/100, f, str);
    t = toc;
end
end

function [frequency] = freq(x, x_vals, y_vals)
n = length(x_vals);
int = x_vals(end) - x_vals(1);
delta = int/(n - 1);

x1 = floor((x - x_vals(1))/delta + x_vals(1));
if x1 <= 0
    x1 = 1;
end
x2 = ceil((x - x_vals(1))/delta + x_vals(1));
if x2 <= 0
    x2 = 1;
end

if x_vals(x2) - x_vals(x1) <= 0
    frequency = y_vals(x1);
else
    frequency = (y_vals(x2) - y_vals(x1))/(x_vals(x2) - x_vals(x1));
    frequency = frequency*(x - x_vals(x1)) + y_vals(x1);
end
end


function [avg_p, mode_p] = trapezoid_data(average, mode, x_vals, y_vals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -Error Checks- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nargin < 4) % Check if less than two inputs
    error(['At least four arguments are required:'...
        'average, mode, x_vals, and y_vals.']);
end

n = length(x_vals); m = length(y_vals);
if(n ~= m) % Check if different # of data
    error('The number of y values must be the same as x values.'); % points
end

%%%%%%%%%%%%%%%%%%%%%%%% -Begin Trapezoidal Data- %%%%%%%%%%%%%%%%%%%%%%%%%
int = x_vals(end) - x_vals(1);
delta = int/(n - 1);
total = trapezoid_data_step(x_vals, y_vals);

avg_ind1 = floor((average - x_vals(1))/delta + x_vals(1)) - x_vals(1) + 1;
avg_ind2 = ceil((average - x_vals(1))/delta + x_vals(1)) - x_vals(1) + 1;
mode_ind = floor((mode - x_vals(1))/delta + x_vals(1)) - x_vals(1) + 1;

if avg_ind1 ~= avg_ind2
    y_avg = (y_vals(avg_ind1) - y_vals(avg_ind2))*(x_vals(avg_ind1) -...
        average)/(x_vals(avg_ind1) - x_vals(avg_ind2)) + y_vals(avg_ind2);
else
    y_avg = y_vals(avg_ind1);
end

integral = (y_avg + y_vals(avg_ind1))*(average - x_vals(avg_ind1))/2;
integral = integral + trapezoid_data_step(x_vals(1:avg_ind1),...
    y_vals(1:avg_ind1));
avg_p = integral/total;

integral = trapezoid_data_step(x_vals(1:mode_ind), y_vals(1:mode_ind));
mode_p = integral/total;
end

function [point, ng] = reverse_trapezoid(x_vals, y_vals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -Error Checks- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nargin < 2) % Check if less than two inputs
    error(['At least two arguments are required:'...
        'proportion, x_vals, and y_vals.']);
end

n = length(x_vals); m = length(y_vals);
if(n ~= m) % Check if different # of data
    error('The number of y values must be the same as x values.'); % points
end

%%%%%%%%%%%%%%%%%%%% -Begin Reverse Trapezoidal Data- %%%%%%%%%%%%%%%%%%%%%
total = trapezoid_data_step(x_vals, y_vals);
integral = @(x) trapezoid_data_step(x_vals(x:end), y_vals(x:end));
ng = [0.5 0.99 0.999 0.9999];

i = n;
for j = 4:-1:1
    while (integral(i)/total) < (1 - ng(j))
        i = i - 1;
    end
    
    area = 1 - integral(i+1)/total - ng(j);
    y = @(x) (y_vals(i) - y_vals(i+1))*(x - x_vals(i+1))...
        /(x_vals(i) - x_vals(i+1)) + y_vals(i+1);
    func = @(x) area - (y_vals(i+1) + y(x))*(x_vals(i+1) - x)/total/2;
    ng(j) = fzero(func, [x_vals(i) x_vals(i+1)]);
    i = i + 1;
end

point = ng(1);
ng = ng(2:4);
end

function [integral] = trapezoid_data_step(x_vals, y_vals)
% Initialize values
n = length(x_vals);
num_slices = n-1;
x_min = x_vals(1);
x_max = x_vals(end);
integral = 0;

if x_max == x_min
    integral = 0;
else
    for i = 1:num_slices
        integral = integral + (x_vals(i+1) - x_vals(i))*(y_vals(i) +...
            y_vals(i+1))/2;
    end
end
end


function [runs] = averageRun(varargin)
probs = cell2mat(varargin); sumP = sum(probs);
if sumP <= 0
    runs = 0;
    return
elseif sumP > 1
    probs = probs/sumP;
    sumP = 1;
end

runs = 1;
for i = 1:length(probs)
    probArray = probs; probArray(i) = [];
    run = averageRun(probArray);
    runs = runs + probs(i)*run;
end

runs = runs/sumP;
end


function numOut = addComma(numIn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%% -Begin Add Comma- %%%%%%%%%%%%%%%%%%%%%%%%%%%%
jf=java.text.DecimalFormat; % comma for thousands, three decimal places
numOut= char(jf.format(numIn)); % omit "char" if you want a string out
end
