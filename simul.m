function simul(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINNER
% July 13th, 2020
% Drop Simulation

% For a specified number of iterations and desired probibilities simul()
% will count the number of times it takes to succeed for each probability
% at least once, plot the data, and calculated the average/expected value
% as well as the values at which 99%, 99.9%, and 99.99% certainty is
% reached.

% Inputs:
% Any list of real, positive probabilities or integers

% Outputs:
% None

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -Error Checks- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nargin % Check for any negative or zero valued chances
    if(varargin{i} <= 0)
        error('Chances must be greater than zero.');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%% -Begin Simulation & Plot- %%%%%%%%%%%%%%%%%%%%%%%%
samples = 1000000; % number of iterations to simulate

n = nargin;
probs = cell2mat(varargin);
if sum(probs) >= 1 % if total probability exeeds 100%, use their weights
    probs = probs / sum(probs);
end
x_vals = 0:n; sample_data = 0; null_prob = 1-sum(probs);

if  null_prob < 0
    null_prob = 0;
end

tic; t = 0; time = 0; del = 1/samples;
str = sprintf(['Progress: 0%%\n',...
    'ETA: 0 Seconds, 0 Minutes, 0 Hours, 0 Days']);
f = waitbar(0, str, 'Name', 'Simulating...', 'CreateCancelBtn',...
    'setappdata(gcbf, ''canceling'', 1)');
setappdata(f, 'canceling', 0);
for k = 1:samples
    if getappdata(f, 'canceling')
        delete(f);
        return
    end
    count = zeros(n + 1, 1);
    while sum(count) < samples
        if getappdata(f, 'canceling')
            delete(f);
            return
        end
        j = randsample(x_vals, 1, true, [null_prob probs]);
        if j >= 1
            count(j) = count(j) + 1;
        else
            count(n + 1) = count(n + 1) + 1;
        end
        
        if sum(count(1:n) >= 1) == n
            break;
        end
    end
    if length(sample_data) < sum(count)
        sample_data((length(sample_data) + 1):sum(count)) =...
            zeros(sum(count)-length(sample_data), 1);
    end
    sample_data(sum(count)) = sample_data(sum(count)) + 1;
    
    if mod(k/samples, del) < 1/samples
        if getappdata(f, 'canceling')
            delete(f);
            return
        end
        
        i = 100*k/samples;
        t_old = t; t = toc; time = time + (t - t_old);
        eta = time/k; eta = eta*(samples - k);
        eta = eta/60/60/24; days = fix(eta); hrs = fix((eta - days)*24);
        mins = fix(((eta - days)*24 - hrs)*60);
        secs = fix((((eta - days)*24 - hrs)*60 - mins)*60);
        eta = [num2str(secs) ' Seconds, ' num2str(mins) ' Minutes, '...
            num2str(hrs) ' Hours, ' num2str(days) ' Days'];
        str = sprintf(['Progress: ', num2str(floor(i)), '%%\n',...
            'ETA: ', eta]);
        waitbar(i/100, f, str);
        
        if (t - t_old) < 0.5
            del = del*2;
        elseif (t - t_old) > 1.5
            del = del/2;
        end
    end
end
waitbar(1, f, 'Simulation Finished!'); pause(1); delete(f);

max_n = length(sample_data);
x_vals = 1:max_n; y_vals = sample_data;
d_x_vals = x_vals(y_vals > 0); d_y_vals = y_vals(y_vals > 0);
average = sum(d_x_vals.*d_y_vals)/samples;

x_vals = linspace(d_x_vals(1), d_x_vals(end), 10000);
y_vals = x_vals;

[~, i] = max(x_vals >= average);
y_vals1 = interp1(d_x_vals, d_y_vals, x_vals, 'spline');
y_vals2 = interp1(d_x_vals, d_y_vals, x_vals, 'pchip');
y_vals(1:i) = y_vals1(1:i); y_vals(i+1:end) = y_vals2(i+1:end);

mode = x_vals(y_vals == max(y_vals));
[average_p, mode_p] = trapezoid_data(average, mode, x_vals, y_vals);

[median, ng] = reverse_trapezoid(d_x_vals, d_y_vals);
range = round((ng(3) - ng(1))/2, 0);

varience = sum(d_y_vals.*(d_x_vals.^2))/samples - average^2;
std_err = sqrt(varience);
se_int = [-std_err; std_err];
average_int(1:2) = 1.96*se_int + average;

time = toc/60/60/24; days = fix(time); hrs = fix((time - days)*24);
mins = fix(((time - days)*24 - hrs)*60);
secs = (((time - days)*24 - hrs)*60 - mins)*60;

%
figure; % from here on is mostly just plotting customization
cmap = colormap(jet(length(x_vals)));
set(gcf, 'position', [10 50 800 600]);
c = linspace(0, 1, length(x_vals)); y_vals = y_vals/samples;
scatter(x_vals, y_vals, 0.5, c); hold on
c = interp1(x_vals, c, d_x_vals); d_y_vals = d_y_vals/samples;
scatter(d_x_vals, d_y_vals, 36, c);
yt = get(gca, 'YTick'); delta = (yt(end) - yt(1))/(length(yt) - 1);

yLim = [0 yt(end)+2*delta]; set(gca, 'YLim', yLim);
str = ['  CI_{50%%}= ' char(177), num2str(0.6745*std_err,'%.5f')...
    '    CI_{75%%}= ' char(177), num2str(1.15*std_err,'%.5f')...
    '    CI_{90%%}= ' char(177), num2str(1.645*std_err,'%.5f')...
    '    CI_{95%%}= ' char(177), num2str(1.96*std_err,'%.5f')...
    '    CI_{99%%}= ' char(177), num2str(2.576*std_err,'%.5f')];
text(n, yLim(2)*0.0275, sprintf(str));

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
xlim([0 max_n]); ax = gca; ax.YGrid = 'on'; ax.XGrid = 'on';
set(gca, 'xscale', 'log');
hold off
%}

fprintf(['   Time Elapsed      = ' num2str(secs) ' Seconds, '...
    num2str(mins) ' Minutes, ' num2str(hrs) ' Hours, ' num2str(days)...
    ' Days\n']);
fprintf(['   Expected          = ' num2str(average_int(1)) ' - '...
    num2str(average_int(2)) ' (' num2str(average) ' - '...
    num2str(average_p*100) '%%)\n']);
fprintf(['   Nearly Guaranteed = ' num2str(round(ng(2))) ' '...
    char(177) ' ' num2str(range) '\n\n']);
end


function [frequency] = freq(x, x_vals, y_vals)
n = length(x_vals);
int = x_vals(end) - x_vals(1);
delta = int/(n - 1);

x1 = floor((x - x_vals(1))/delta + x_vals(1));
x2 = ceil((x - x_vals(1))/delta + x_vals(1));

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

avg1 = floor((average - x_vals(1))/delta + x_vals(1));
avg2 = ceil((average - x_vals(1))/delta + x_vals(1));
mode1 = floor((mode - x_vals(1))/delta + x_vals(1));

y_avg = (y_vals(avg2) - y_vals(avg1))*(average - x_vals(avg1))...
    /(x_vals(avg2) - x_vals(avg1)) + y_vals(avg1);
integral = (y_avg + y_vals(avg1))*(average - x_vals(avg1))/2;
integral = integral + trapezoid_data_step(x_vals(1:avg1), y_vals(1:avg1));
avg_p = integral/total;

integral = trapezoid_data_step(x_vals(1:mode1), y_vals(1:mode1));
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
    while integral(i)/total < 1 - ng(j)
        i = i - 1;
    end
    area = 1 - ng(j) - integral(i+1)/total;
    y = @(x) (y_vals(i+1) - y_vals(i))*(x_vals(i+1) - x)...
        /(x_vals(i+1) - x_vals(i)) + y_vals(i+1);
    func = @(x) (y_vals(i+1) + y(x))*(x_vals(i+1) - x)/total - area;
    ng(j) = IQI(func, x_vals(i), (x_vals(i)+x_vals(i+1))/2, x_vals(i+1));
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


function [x_root,func_val,error_approx,num_iterations] = IQI(func,...
    x_guess_one,x_guess_two,x_guess_three,error_desired,max_iterations)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -Error Checks- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nargin < 4) % Check if less than four inputs
    error('At least four arguments are required.');
end
if(x_guess_one == x_guess_two) % Check if both guesses are the same
    error('x_guess_one must not be equal to x_guess_two.');
end
if(x_guess_one == x_guess_three) % Check if both guesses are the same
    error('x_guess_one must not be equal to x_guess_three.');
end
if(x_guess_two == x_guess_three) % Check if both guesses are the same
    error('x_guess_two must not be equal to x_guess_three.');
end
if(nargin < 5) || isempty(error_desired) % Default error_desired if not
    error_desired = 0.0001;              % given
end
if(nargin < 6) || isempty(max_iterations) % Default max_iterations if not
    max_iterations = 50;                  % given
end
if(max_iterations <= 0) % Check if max_iterations is greater than zero
    error('max_iterations must be greater than zero.');
end

%%%%%%%%%%%%%%%%% -Begin Inverse Quadratic Interpolation- %%%%%%%%%%%%%%%%%
% Initialize values
x_old_old = x_guess_one; x_old = x_guess_two; x_current = x_guess_three;
error_approx = error_desired + 1;
num_iterations = 0;

if(func(x_old_old) == 0)
    x_root = x_old_old;
    func_val = func(x_old_old);
    error_approx = 0;
    num_iterations = 0;
elseif(func(x_old) == 0)
    x_root = x_old;
    func_val = func(x_old);
    error_approx = 0;
    num_iterations = 0;
elseif(func(x_current) == 0)
    x_root = x_current;
    func_val = func(x_current);
    error_approx = 0;
    num_iterations = 0;
else
    while(error_approx > error_desired)&& (num_iterations < max_iterations)
        num_iterations = num_iterations + 1;
        if(func(x_current) ~= func(x_old)&& func(x_old) ~=...
                func(x_old_old)&& func(x_current) ~= func(x_old_old))
            x_new = (func(x_old)*func(x_old_old)/((func(x_current)-...
                func(x_old))*(func(x_current)-func(x_old_old))))*...
                x_current+(func(x_current)*func(x_old_old)/((func(x_old)...
                -func(x_current))*(func(x_old)-func(x_old_old))))*x_old+...
                (func(x_current)*func(x_old)/((func(x_old_old)-...
                func(x_current))*(func(x_old_old) - func(x_old))))*...
                x_old_old;
        else
            if(func(x_current) == func(x_old))
                error(['x_current and x_old have equal function values,'...
                    ' resulting in a division by zero. Try a different'...
                    ' guess two or guess three.'])
            elseif(func(x_current) == func(x_old_old))
                error(['x_current and x_old_old have equal function'...
                    ' values, resulting in a division by zero. Try a'...
                    ' different guess one or guess three.'])
            else
                error(['x_old and x_old_old have equal function'...
                    ' values, resulting in a division by zero. Try a'...
                    ' different guess one or guess two.'])
            end
        end
        if(x_new ~= 0)
            error_approx = abs((x_new - x_current)/x_new) * 100;
        end
        x_old_old = x_old;
        x_old = x_current;
        x_current = x_new;
    end
    x_root = x_new;
    func_val = func(x_root);
end

end

function numOut = addComma(numIn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%% -Begin Add Comma- %%%%%%%%%%%%%%%%%%%%%%%%%%%%
jf=java.text.DecimalFormat; % comma for thousands, three decimal places
numOut= char(jf.format(numIn)); % omit "char" if you want a string out

end
