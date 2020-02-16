function simul(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINNER
% February 2, 2020
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
x = 0:n; sum_data = 0; sample_data = 0; null_prob = 1-sum(probs);

if  null_prob < 0
    null_prob = 0;
end

tic
for k = 1:samples
    count = zeros(n + 1, 1);
    while sum(count) < samples
        j = randsample(x, 1, true, [null_prob probs]);
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
            zeros(1, sum(count)-length(sample_data));
    end
    sample_data(sum(count)) = sample_data(sum(count)) + 1;
end
for i=1:length(sample_data)
    if length(sum_data) < i
        sum_data((length(sum_data)+1):i) = zeros(1, i-length(sum_data));
    end
    sum_data(i) = sum_data(i) + sample_data(i);
end

max_n = length(sum_data);
x_vals = 1:max_n; y_vals = sum_data;

average = sum(x_vals.*y_vals)/samples;
average_p = trapezoid_data(average, x_vals, y_vals, n);
ng(1) = reverse_trapezoid(0.99, x_vals, y_vals, n);
ng(2) = reverse_trapezoid(0.999, x_vals, y_vals, n);
ng(3) = reverse_trapezoid(0.9999, x_vals, y_vals, n);
range = round((ng(3) - ng(1))/2,0);

varience = sum((x_vals.*y_vals - average).^2)/(samples -  1);
std_err = sqrt(varience/samples);
se_int = [-std_err; std_err];
average_int(1:2) = 2*se_int + average;

time = toc/60/60/24; days = fix(time); hrs = fix((time - days)*24);
mins = fix(((time - days)*24 - hrs)*60);
secs = (((time - days)*24 - hrs)*60 - mins)*60;

% from here on is mostly just plotting customization
x_vals = x_vals(y_vals > 0); y_vals = y_vals(y_vals > 0)/samples;
c = linspace(0, 1, length(x_vals));

figure;
cmap = colormap(cool(max_n)); set(gcf, 'position', [10 50 800 600]);
scatter(x_vals, y_vals, 36, c); hold on
yLim = get(gca, 'YLim');

plot([average average], yLim, 'Color', cmap(round(average), :));
text(average + max_n*0.01, max(y_vals)*0.85, sprintf(['Average\n'...
    num2str(average,'%.4f') ' ' char(177) ' ' num2str(2*std_err,'%.5f')...
    '\n(' num2str(average_p*100,'%.3f') '%%)']));

plot([ng(1) ng(1)], yLim, 'Color', cmap(round(ng(1)), :));
text(ng(1) + max_n*0.01, max(y_vals)*0.85, sprintf([num2str(ng(1),'%.3f')...
    '\n(99%%)']));

plot([ng(2) ng(2)], yLim, 'Color', cmap(round(ng(2)), :));
text(ng(2) + max_n*0.01, max(y_vals)*0.85, sprintf([num2str(ng(2),'%.3f')...
    '\n(99.9%%)']));

plot([ng(3) ng(3)], yLim, 'Color', cmap(round(ng(3)), :));
text(ng(3) + max_n*0.01, max(y_vals)*0.85, sprintf([num2str(ng(3),'%.3f')...
    '\n(99.99%%)']));

title(['Number of Runs for ' addComma(samples) ' Samples']);
xlabel('Number of Runs'); ylabel('Frequency');
xlim([0 max_n]); ax = gca; ax.YGrid = 'on'; ax.XGrid = 'on';
hold off

fprintf(['   Time Elapsed      = ' num2str(secs) ' Seconds, '...
    num2str(mins) ' Minutes, ' num2str(hrs) ' Hours, ' num2str(days)...
    ' Days\n']);
fprintf(['   Expected          = ' num2str(average_int(1)) ' - '...
    num2str(average_int(2)) ' (' num2str(average) ' - '...
    num2str(average_p*100) '%%)\n']);
fprintf(['   Nearly Guaranteed = ' num2str(round(ng(2))) ' '...
    char(177) ' ' num2str(range) '\n\n']);
end


function [proportion] = trapezoid_data(point, x_vals, y_vals, n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -Error Checks- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nargin < 4) % Check if less than two inputs
    error(['At least four arguments are required,'...
        'desired, x_vals, y_vals, n']);
end

[x_rows, x_cols] = size(x_vals); [y_rows, y_cols] = size(y_vals);
if(x_rows ~= y_rows || x_cols ~= y_cols) % Check if different # of data
    error('The number of y values must be the same as x values.'); % points
end

%%%%%%%%%%%%%%%%%%%%%%%% -Begin Trapezoidal Data- %%%%%%%%%%%%%%%%%%%%%%%%%
total = trapezoid_data_step(x_vals, y_vals, n);
x1 = floor(point); x2 = ceil(point);
y_point = (y_vals(x2) - y_vals(x1))*(point - x1) + y_vals(x1);
integral = (y_point + y_vals(x1))*(point - x1)/2;
integral = integral + trapezoid_data_step(x_vals(1:x1), y_vals(1:x1), n);
proportion = integral/total;

end

function [point] = reverse_trapezoid(proportion, x_vals, y_vals, n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -Error Checks- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nargin < 4) % Check if less than two inputs
    error(['At least four arguments are required,'...
        'desired, x_vals, y_vals, n']);
end

[x_rows, x_cols] = size(x_vals); [y_rows, y_cols] = size(y_vals);
if(x_rows ~= y_rows || x_cols ~= y_cols) % Check if different # of data
    error('The number of y values must be the same as x values.'); % points
end

%%%%%%%%%%%%%%%%%%%% -Begin Reverse Trapezoidal Data- %%%%%%%%%%%%%%%%%%%%%
total = trapezoid_data_step(x_vals, y_vals, n);
integral = @(x) trapezoid_data_step(x_vals(1:x), y_vals(1:x), n);

for i = n:length(x_vals)
    if integral(i)/total >= proportion
        area = proportion - integral(i-1)/total;
        y = @(x) (y_vals(i-1) - y_vals(i))*x + y_vals(i-1);
        func = @(x) (y_vals(i-1) + y(x-i+1))*(x-i+1)/(2*total) - area;
        point = IQI(func, i-1, i-0.5, i);
        break
    end
end

end

function [integral] = trapezoid_data_step(x_vals, y_vals, n)
% Initialize values
rows = length(x_vals);
num_slices = rows-1;
x_min = x_vals(n);
x_max = x_vals(end);
integral = 0;

if x_max == x_min
    integral = 0;
else
    for i = n:num_slices
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
