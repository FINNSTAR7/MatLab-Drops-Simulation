function simul(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINNER
% March 30, 2019
% Drop Simulation

% For a specified number of iterations and desired probibilities simul()
% will count the number of times it takes to succeed for each probability
% at least once, plot the data, and calculate the average/expected value
% as well as the values at which 99%, 99.9%, and 99.99% certainty is
% reached.

% Inputs:
% Any numeric list of real, positive probabilities or integers

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

n = nargin; probs(1:n,1) = cell2mat(varargin);
if sum(probs) > 1 % if total probability exeeds 100%, use their weights
    probs = probs / sum(probs);
end
x = 0:n; sum_data = 0;

tic
sample_data = 0;
for k=1:samples
    count = zeros(n + 1,1);
    while sum(count) < samples
        b = 0;
        j = randsample(x, 1, true, [1-sum(probs) probs']);
        if j >= 1
            count(j) = count(j) + 1;
        else
            count(n + 1) = count(n + 1) + 1;
        end
        
        for i=1:n
            if count(i) >= 1
                b = b + 1;
            end
        end
        if b == n
            break;
        end
    end
    if length(sample_data) < sum(count)
        sample_data((length(sample_data)+1):sum(count)) =...
            zeros(1,sum(count)-length(sample_data));
    end
    sample_data(sum(count)) = sample_data(sum(count)) + 1;
end
for i=1:length(sample_data)
    if length(sum_data) < i
        sum_data((length(sum_data)+1):i) = zeros(1,i-length(sum_data));
    end
    sum_data(i) = sum_data(i) + sample_data(i);
end

i = length(sum_data); avg_data = zeros(i,2);
avg_data(1:i,1) = 1:i; avg_data(1:i,2) = sum_data;

data = avg_data; max = size(data,1);
x_vals = data(1:max,1); y_vals = data(1:max,2);

average = sum(prod(data(:,1:2),2)) / sum(y_vals);
average_p = trapezoid_data(average,x_vals,y_vals,n);

time = toc/60/60/24; days = fix(time); hrs = fix((time - days)*24);
mins = fix(((time - days)*24 - hrs)*60);
secs = (((time - days)*24 - hrs)*60 - mins)*60;

% from here on is mostly just plotting customization
sz = zeros(x_vals(max)); c = linspace(0,1,x_vals(max));
for i=1:x_vals(max)
    if y_vals(i) <= 0
        sz(i) = 1;
        c(i) = 0;
    else
        sz(i) = 36;
    end
end
mag = [1, 0, 1]; cyan = [0, 1, 1]; figure
colors_p = [linspace(mag(1),cyan(1),max)',...
    linspace(mag(2),cyan(2),max)',...
    linspace(mag(3),cyan(3),max)'];
cmap = colormap(colors_p); set(gcf,'position',[10,50,800,600]);
scatter(x_vals,y_vals/sum(y_vals),sz(:,1),c); hold on

vline(average,cmap(round(average),:),sprintf(['Average\n'...
    num2str(average_p*100,'%.3f') '%%\n(' num2str(average,'%.4f')...
    ')']),0.92,0.01);

ng(1,1) = reverse_trapezoid(0.99,x_vals,y_vals,n);
ng(1,2) = reverse_trapezoid(0.999,x_vals,y_vals,n);
ng(1,3) = reverse_trapezoid(0.9999,x_vals,y_vals,n);

vline(ng(1,1),cmap(round(ng(1,1)),:),sprintf(['99%%\n('...
    num2str(ng(1,1),'%.3f') ')']),0.75,0.01);
vline(ng(1,2),cmap(round(ng(1,2)),:),sprintf(['99.9%%\n('...
    num2str(ng(1,2),'%.3f') ')']),0.9,0.01);
vline(ng(1,3),cmap(round(ng(1,3)),:),sprintf(['99.99%%\n('...
    num2str(ng(1,3),'%.3f') ')']),0.75,0.01);

title(['Number of Runs for ' addComma(sum(y_vals)) ' Samples']);
xlabel('Number of Runs'); ylabel('Frequency');
int = linspace(n,max,ceil(max/n));
xlim([0 max]); xticks(ceil(int)); ax = gca; ax.YGrid = 'on';
ax.XGrid = 'on';
hold off

range = round((ng(1,3) - ng(1,1))/2,0);
fprintf(['   Time Elapsed      = ' num2str(secs) ' Seconds, '...
    num2str(mins) ' Minutes, ' num2str(hrs) ' Hours, ' num2str(days)...
    ' Days\n']);
fprintf(['   Expected          = ' num2str(floor(average)) ' - '...
    num2str(ceil(average)) ' (' num2str(average) ' - '...
    num2str(average_p*100) '%%)\n']);
fprintf(['   Nearly Guaranteed = ' num2str(round(ng(1,2),0)) ' '...
    char(177) ' ' num2str(range) '\n\n']);
end


function [integral] = quad_adapt(func,x_min,x_max,tolerance)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -Error Checks- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nargin < 3) % Check if less than three inputs
    error(['At least three arguments are required, a function, x min,'...
        ' and x max.']);
end
if(x_max <= x_min) % Check if min is greater than max, or are the same
    error('x_max must be greater than x_min.');
end
if(nargin < 4) || isempty(tolerance) % Default tolerance if not
    tolerance = 0.000001;            % given
end

%%%%%%%%%%%%%%%%%%%%%%% -Begin Adaptive Quadrature- %%%%%%%%%%%%%%%%%%%%%%%
integral = quad_step(func,x_min,x_max,tolerance);

end

function [integral] = quad_step(func,x_min,x_max,tolerance)
% Initialize values
x_bar = (x_min+x_max)/2;
x_barL = (x_min+x_bar)/2; x_barR = (x_bar+x_max)/2;
h_1 = (x_max-x_min)/2; h_2 = (x_max-x_min)/4;

int_1 = h_1*(func(x_min)+4*func(x_bar)+func(x_max))/3;
int_2 = h_2*(func(x_min)+4*func(x_barL)+2*func(x_bar)+4*func(x_barR)+...
    func(x_max))/3;

if abs(int_1 - int_2) <= tolerance
    integral = (16*int_2 - int_1)/15;
else
    intL = quad_step(func,x_min,x_bar,tolerance);
    intR = quad_step(func,x_bar,x_max,tolerance);
    integral = intL + intR;
end

end


function [proportion] = trapezoid_data(point,x_vals,y_vals,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -Error Checks- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nargin < 4) % Check if less than two inputs
    error(['At least four arguments are required,'...
        'desired, x_vals, y_vals, n']);
end

[x_rows,x_cols] = size(x_vals); [y_rows,y_cols] = size(y_vals);
if(x_rows ~= y_rows || x_cols ~= y_cols) % Check if different # of data
    error('The number of y values must be the same as x values.'); % points
end

%%%%%%%%%%%%%%%%%%%%%%%% -Begin Trapezoidal Data- %%%%%%%%%%%%%%%%%%%%%%%%%
total = trapezoid_data_step(x_vals,y_vals,n);
x1 = floor(point); x2 = ceil(point);
y_point = (y_vals(x2) - y_vals(x1))*(point - x1) + y_vals(x1);
integral = (y_vals(x1) + y_point)*(point - x1)/2;
integral = integral + trapezoid_data_step(x_vals(1:x1,1),y_vals(1:x1,1),n);
proportion = integral/total;

end

function [point] = reverse_trapezoid(proportion,x_vals,y_vals,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -Error Checks- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nargin < 4) % Check if less than two inputs
    error(['At least four arguments are required,'...
        'desired, x_vals, y_vals, n']);
end

[x_rows,x_cols] = size(x_vals); [y_rows,y_cols] = size(y_vals);
if(x_rows ~= y_rows || x_cols ~= y_cols) % Check if different # of data
    error('The number of y values must be the same as x values.'); % points
end

%%%%%%%%%%%%%%%%%%%% -Begin Reverse Trapezoidal Data- %%%%%%%%%%%%%%%%%%%%%
total = trapezoid_data_step(x_vals,y_vals,n);
integral = @(x) trapezoid_data_step(x_vals(1:x,1),y_vals(1:x,1),n);
for i=n:x_rows
    if integral(i)/total >= proportion
        area = proportion - integral(i-1)/total;
        y = @(x) (y_vals(i-1,1)-y_vals(i,1))*x + y_vals(i-1,1);
        func = @(x) (y_vals(i-1,1)+y(x-i+1))*(x-i+1)/(2*total) - area;
        point = IQI(func,i-1,i-0.5,i);
        break
    end
end

end

function [integral] = trapezoid_data_step(x_vals,y_vals,n)
% Initialize values
column_size = size(x_vals,1);
num_slices = column_size-1;
x_min = x_vals(n,1);
x_max = x_vals(column_size,1);
integral = 0;

if x_max == x_min
    integral = 0;
else
    for i = n:num_slices
        integral = integral+(x_vals(i+1,1)-x_vals(i,1))*(y_vals(i,1)+...
            y_vals(i+1,1))/2;
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


function hhh=vline(x,in1,in2,high,tall)
% function h=vline(x, linetype, label)
%%%%%%%%%%%%%%%%%%%%%%%%%% -Begin Vertical Line- %%%%%%%%%%%%%%%%%%%%%%%%%%
if length(x)>1  % vector input
    for I=1:length(x)
        switch nargin
            case 1
                linetype='r:';
                label='';
            case 2
                if ~iscell(in1)
                    in1={in1};
                end
                if I>length(in1)
                    linetype=in1{end};
                else
                    linetype=in1{I};
                end
                label='';
            case 3
                if ~iscell(in1)
                    in1={in1};
                end
                if ~iscell(in2)
                    in2={in2};
                end
                if I>length(in1)
                    linetype=in1{end};
                else
                    linetype=in1{I};
                end
                if I>length(in2)
                    label=in2{end};
                else
                    label=in2{I};
                end
        end
        h(I)=vline(x(I),linetype,label);
    end
else
    switch nargin
        case 1
            linetype='r:';
            label='';
        case 2
            linetype=in1;
            label='';
        case 3
            linetype=in1;
            label=in2;
        case 4
            linetype=in1;
            label=in2;
        case 5
            linetype=in1;
            label=in2;
    end
    
    g=ishold(gca);
    hold on
    
    y=get(gca,'ylim');
    h=plot([x x],y,'Color',linetype);
    if length(label)
        xx=get(gca,'xlim');
        xrange=xx(2)-xx(1);
        text(x+tall*xrange,y(1)+high*(y(2)-y(1)),label,'color',[0 0 0])
    end
    
    if g==0
        hold off
    end
    set(h,'tag','vline','handlevisibility','off')
end % else

if nargout
    hhh=h;
end

end


function numOut = addComma(numIn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%% -Begin Add Comma- %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    jf=java.text.DecimalFormat; % comma for thousands, three decimal places
    numOut= char(jf.format(numIn)); % omit "char" if you want a string out
    
end
