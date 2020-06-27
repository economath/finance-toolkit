% Portfolio Optimization Using Linear Algebra
% by Vishak Srikanth
% The toolkits uses linear algebra to provide solutions for 3 common scenarios for portfolio optimization
% Given a number of securities with varied returns a1, a2,..., an, investors are interested in determing how they can
% maximize returns while minimizng risk by constructing a blend of securities that can provide them good returns at
% acceptable risk. The Efficient frontier describes the maximum possible expected return for any given amount of risk
% from the portfolio of assets Or alternatively, the minimum amount of risk that one must accept for
% any given amount of expected return. This program helps users construct different portfolios:
% 1. Minimum variance portfolio: Since there is no rewards without risk, rational investors want to mimize risk
% with their security investements. The portfolio variance represents this risk, they want to create a portfolio
% constructed with n securities that provides the minimum global variance.  As you change the weights of the different
% the individual securities, the variance of the overall protfolio changes and investors want to understand what is the
% portfolio that provides the minimum global variance across the spectrum of returns and mixes of given assets.
% 2. Tangency or Risk free asset mix: The tangency portfolio is the portfolio with the highest Sharpe ratio
% shows how much return you can expect to get for every unit of risk. Tangency portfolio represents the market return.
% It’s called the tangency because it’s located at the tangency point of the Capital Market Line and the
% Efficient Frontier.
% 3. Target portfolio with a given return that lies on the efficient frontier.
% This toolkit can create the optimal portfolios mentioned above for any number of assets!
% To create the portfolios the only input data needed is a csv (comma separated values) file where the columns
% are the securities and the rows are the adjusted closing prices on specific dates
% on each date which can be easily obtained from sources like Yahoo Finance
% For testing I used a website provided by Jason Strimpel http://finance.jasonstrimpel.com/bulk-stock-download/
% where you can bulk download for free the pricing data for multiple stocks over the past 60 or so years!
% The adjusted close prices include any dividend/ split information on a specific date
% As long as you have prices for any time interval for all assets we can compute the optimal portfolio!
% This program has been tested with portfolios as large as all S&P 500 stocks to a handful and for daily, weekly,
% monthly data over short and long durations (1 year, 5 year, 10 year, 20 year and 30 year time horizons.
% The best example is a portfolio made only from the DOW JONES 30 stocks monthly reurns over 10 years which illustrate
% how the program works
% Note: Most csv files from the sources such as Yahoo Finance have 1st column as date which can be ignored
% in the calculations as we are only focused on calculating returns across each interval (day, week or month, typically)
% The securities individual returns coputed between two dates is used to determine means, variances, covariances etc.
% that determine the optimal portfolios.


% Need to load various packages for reading csv files and using dataframes and tables
pkg load io
pkg load dataframe

clear
clc

% Users can choose to upload a csv file of adjusted close price data
[file,path,indx] = uigetfile('*.csv',...
               'Select an icon file','icon.png')
if isequal(file,0)
   disp('User selected Cancel')
else
   disp(['User selected ', fullfile(path, file),... 
         ' and filter index: ', num2str(indx)])
end

% Read the csv data file
stock_data = csvread(fullfile(path,file), 1,1);

%stock_df = dataframe(fullfile(path,file))


% Check the data size to see how many stocks and how many data points.
size(stock_data)
%c = dataframe(fullfile(path,file));

% Read the data into a cell array
c= csv2cell(fullfile(path,file));

% c = readtable(fullfile(path,file)','ReadRowNames',true);
%c(1:5,1:4)
%c(1,:) % will print all headers
%c(2:end,:);  % will all  data rows

% Since the 1st column is the date, we start with columns 2-end which hold the data for the stocks
% For example, the DOW JONES DATA SET has the current list of 30 stocks in DJIA
% AAPL	AXP	BA	CAT	CSCO CVX DIS DOW GS	HD IBM INTC JNJ	JPM	KO MCD	MMM	MRK	MSFT NKE PFE PG	RTX	TRV	UNH	V VZ WBA WMT XOM
stock_symb = c(1, 2:end);

size(stock_symb); % get number of stocks needed in portfolio
%dates =c(:, 1)
%stock_data_clean = c(2:end, 2:end);
%stock_data_clean  =  fillmissing(stock_data, "spline")
%stock_data_clean  =  rmmissing(stock_data, 2);

% Since matrix inversion algortihms have issues when we have missing data leading to sparse
% matrices that might be singular, conservatively we can drop all stocks which are missing date for the entire duration
% For example if we take the Dow Jones 30 stock data from 2010-2020 only 1 stock Dow Chemicals has missing data and will
% be dropped. Even with S&P 500 portfolios, over 450 stocks are retained.
% the 'all' command looks for columns where all cells are non-empty
cols_to_keep = find(all(stock_data));

% Create a clean data set after dropping columns or stocks with missing data during the period.
stock_data_clean = stock_data(:, cols_to_keep);
stock_symb_clean = stock_symb(cols_to_keep);

% Approach: We first need to compute the return over each period at tn and tn+1
% return(n) = (Pn+1-Pn)/Pn where Pn = adj. closing price at tn etc.
% For this we use the nifty diff function that does a row shift and computes the finite difference!
diff_mat = -1 .* diff(stock_data_clean);
% Then we measure the return on the basis of the previous period's closing price
stock_returns = diff_mat ./ stock_data_clean(2:end, :);

% Compute the mean, std. dev and covariance matrix for all the stocks
% Note: the diagonal entries of covariance matrix will have the stock's standard deviation
mean_returns = mean(stock_returns);
stdev_returns =  std(stock_returns);
covar_returns =cov(stock_returns, 1);
[covar_nrows, covar_ncols] = size(covar_returns)
[mean_ret_nrows, mean_ret_ncols] = size(mean_returns)
[stdev_ret_nrows, stdev_ret_ncols] = size(stdev_returns)
nstocks = mean_ret_ncols;

% Now comes portfolio computations as outlined in the book Modern Portfolio Theory and Investment
% Analysis by Elton, E., M. Gruber, S. Brown, and W. Goetzmann (2009) - John Wiley & Sons
% The variance-covariance matrix of the assets and their expected returns are key inputs to design the portfolio
% The variance-covariance matrix (Σ) is the symmetric matrix sigma below and the expected or mean returns of
% the individual securtities calculated from the average over the time period (daily, weekly, or monthly for example)
% from above are stored in a column-vector for easier manipulation (zbar).
% We need to calculate the weights that will yield a portfolio with minimum risk.
% The risk of the portfolio is related to the standard deviation and it can be calculated as follows
% So we now have to minimize:
% portfolio variance σ = √X'TΣX subject to: ∑iXi =1, where X is a vector of weights in portfolio
% The expected return of the portfolio is given by w = X'(zbar)
% The weights Xi are fractions invested in each security and thus they must add up to 1.
% The minimization problem can be solved using Lagrange multipliers to minimize the variance.
sigma = covar_returns;
% a = is just a 1's column vector
a=ones(nstocks,1);
% We create a singleton zero
zero_scal = zeros(1);
last_row = [a'  0];
%
A=[2*sigma  a; a'  0];
b=eye(nstocks+1)(:, nstocks+1);
x = A\b;

% Checked multiple ways of computing inverse and premultiplying , best seems to be just the \ operator
%y = inv(A)*b;
%d = det(A);
%x =linsolve(A,b)
%symb_col = string(stock_symb_clean)'
%pf = dataframe(symb_col)
%pf(:, 2) =  x(1:end-1)






% Mean Variance Optimizer
% Next we look at how to compute the efficient frontier and the tangency portfolio
% S is sigma matrix of security variance - covariance of returns
S = covar_returns;
% 
% Vector of security expected returns
zbar = mean_returns';

%Risk Free Asset Return: use latest 3-month T-bills rate 0.125%
Rf = 0.00125

% 
% Unity vector..must have same length as zbar
unity = ones(length(zbar),1);
% 
% The standard deviations of the security are simply diagonals of variance-covariance matrix
stdevs = sqrt(diag(S));
% 
%% Calculate Efficient Frontier
% This calculations below explained in detail in Prof. Eric Zivot's Portfolio Theory Class notes here
% https://faculty.washington.edu/ezivot/econ424/portfolioTheoryMatrix.pdf

A = unity'*S^-1*unity;
B = unity'*S^-1*zbar;
C = zbar'*S^-1*zbar;
D = A*C-B^2;
% 
%% Efficient Frontier
npoints = 100;
mu = linspace(min(zbar),max(zbar),npoints) ;
% 
%% Plot Efficient Frontier
minvar = ((A*mu.^2)-2*B*mu+C)/D;
minstd = sqrt(minvar);
% 
mu_max = max(mu)
mu_min = min(mu)
minstd_max= max(minstd)
minstd_min = min(minstd)

% Let's check the optimal portfolio for a target Return of 15%
mu_tar = 0.15

% Calculate Lambda and Gamma
lambda_target = (C - mu_tar*B)/D;
gamma_target =  (mu_tar*A-B)/D;

loglog(minstd.*100,mu.*100,stdevs.*100,zbar.*100,'*')
title("Efficient Frontier with Portolio of Individual Securities",'fontsize',18)
ylabel('Expected Return (%)','fontsize',18)
xlabel('Standard Deviation (%)','fontsize',18)
set(gca,'xlim',[0, 50]) ;
set(gca,'ylim',[0, 25]) ;
grid ;
%
%

% Mean and Variance of Global Minimum Variance Portfolio
% subscript g for global portfolio values
mu_g = B/A
var_g = 1/A
std_g = sqrt(var_g)
 
% Minimum Variance Portfolio Weights
w_g = (S^-1*unity)/A
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tangency Portfolio with a Risk Free Asset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Expected Return of Tangency Portfolio
ztan = (C-B*Rf)/(B-A*Rf);

% Variance and Standard Deviation of Tangency Portfolio
vartan = (C-2*Rf*B + Rf^2*A)/((B-A*Rf)^2);
stdtan = sqrt(vartan);

% Weights for Tangency Portfolio
w_tan = (S^-1*(zbar - Rf*unity))/(B-A*Rf)

% Tangency Line
mu_tan = mu(mu >= Rf);
minvar_rf = (mu_tan-Rf).^2/(C-2*Rf*B+A*Rf^2);
minstd_rf = sqrt(minvar_rf);

% Weights for w_d (tangency when R=0)
w_d = (S^-1*zbar)/B;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Portfolio with Expected Return = 15%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Weights for portfolio with target return = 15%
 
w_s = (lambda_target*A)*w_g + (gamma_target*B)*w_d;
 
% Expected Return of Target Portfolio (should match target)
mu_s = w_s'*zbar
 
% Variance and Standard Deviation of target portfolio
var_s = w_s'*S*w_s;
std_s = sqrt(var_s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Target Return Portfolio w/and w/o Risk Free Asset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subscript s is for target portfolio with an expected return chosen by user
% With and without the presence of the risk free asset the porfolio returns and stdevs are calculated
% Weights for portfolio with target return = 15%, w/o risk-free asset

w_s = (lambda_target*A)*w_g + (gamma_target*B)*w_d;

% Expected Return of Target Portfolio (should match target: 15%)
mu_s = w_s'*zbar;

% Variance and Standard Deviation of target portfolio
var_s = w_s'*S*w_s;
std_s = sqrt(var_s);

% Weights for portfolio with target return = 15%, w/risk free asset

y = (mu_tar - Rf)/(ztan-Rf);
stdtar = stdtan*y;
















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Efficient Frontier and Key Portfolio Points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
figure
%loglog(100.*minstd_rf,'linewidth',2,100.*mu_tan,minstd.*100,mu.*100,'linewidth',2,std_d.*100,mu_d.*100,'*','linewidth',2,std_g.*100,mu_g.*100,'x','linewidth',2,std_s.*100,mu_s.*100,'o','linewidth',2)

lgplot = loglog(100.*minstd_rf,'linewidth',2,100.*mu_tan,100.*minstd,'linewidth',2,100.*mu,100.*stdtan,100.*ztan,'*','linewidth',2,100.*std_g,100.*mu_g,'x','linewidth',2,100.*std_s,100.*mu_s,'x','linewidth',2,100.*stdtar,100.*mu_tar,'*','linewidth',2)
text(.1,100.*Rf,'RF','fontsize',12);
text(0.5+std_g.*100,mu_g.*100,'Global Minimum Variance Portfolio','fontsize',14);
%text(0.5+std_d.*100,mu_d.*100,'Tangency Portfolio when R=0','fontsize',14);
text(0.5+100*stdtan,100*ztan,'Tangency Portfolio','fontsize',12);
text(0.5+std_s.*100,mu_s.*100,'Portfolio with Target Return of 15%','fontsize',14);
%text(100.*stdtar-8,mu_tar.*100+0.5,'Target Return of 15% w/Risk-Free Asset','fontsize',12);
title('Efficient Frontier with Selected Portfolios','fontsize',18)
ylabel('Expected Return (%)','fontsize',18)
xlabel('Standard Deviation (%)','fontsize',18)
waitfor(lgplot);





