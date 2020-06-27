## finance-toolkit
This project is a repository for finance toolkit using Linear Algebra and Octave that I created for Linear Algebra class

 # Portfolio Optimization Using Linear Algebra
 ### by Vishak Srikanth
 
 This toolkit uses Linear Algebra to provide solutions for 3 common scenarios for portfolio optimization
 
 Given a number of securities with varied returns a_1, a_2,..., a_n investors are interested in determing how they can maximize returns while minimizng risk by constructing a blend of securities that can provide them good returns at acceptable risk. 
 
 The Efficient frontier describes the maximum possible expected return for any given amount of risk from the portfolio of assets Or alternatively, the minimum amount of risk that one must accept for any given amount of expected return. 
 
 This program helps users construct different portfolios:
 
 1. Minimum variance portfolio: Since there is no rewards without risk, rational investors want to mimize risk with their security investements. The portfolio variance represents this risk, they want to create a portfolio constructed with n securities that provides the minimum global variance.  As you change the weights of the different the individual securities, the variance of the overall protfolio changes and investors want to understand what is the portfolio that provides the minimum global variance across the spectrum of returns and mixes of given assets.
 
 2. Tangency or Risk free asset mix: The tangency portfolio is the portfolio with the highest Sharpe ratio shows how much return you can expect to get for every unit of risk. Tangency portfolio represents the market return.  It’s called the tangency because it’s located at the tangency point of the Capital Market Line and the Efficient Frontier.
 
 3. Target portfolio with a given return that lies on the efficient frontier. This toolkit can create the optimal portfolios mentioned above for any number of assets!
 
 
## How to create and visualize the portfolios

### Inputs 
The only input data needed is a csv (comma separated values) file where the columns are the securities and the rows are the adjusted closing prices on specific dates on each date which can be easily obtained from sources like Yahoo Finance

## Testing
For testing I used a website provided by Jason Strimpel http://finance.jasonstrimpel.com/bulk-stock-download/ where you can bulk download for free the pricing data for multiple stocks over the past 60 or so years!

The adjusted close prices include any dividend/ split information on a specific date

As long as you have prices for any time interval for all assets we can compute the optimal portfolio!

This program has been tested with portfolios as large as all S&P 500 stocks to a handful and for daily, weekly, monthly data over short and long durations (1 year, 5 year, 10 year, 20 year and 30 year time horizons.

The best example is a portfolio made only from the DOW JONES 30 stocks monthly reurns over 10 years which illustrate how the program works

### Note: 
Most csv files from the sources such as Yahoo Finance have 1st column as date which can be ignored in the calculations as we are only focused on calculating returns across each interval (day, week or month, typically)

The securities individual returns computed between two dates is used to determine means, variances, covariances etc. that determine the optimal portfolios.

## Portfolio Theory

Now comes portfolio computations as outlined in the book Modern Portfolio Theory and Investment Analysis by Elton, E., M. Gruber, S. Brown, and W. Goetzmann (2009- John Wiley & Sons

The variance-covariance matrix of the assets and their expected returns are key inputs to design the portfolio

The variance-covariance matrix (Σ) is the symmetric matrix sigma below and the expected or mean returns of the individual securtities calculated from the average over the time period (daily, weekly, or monthly for example) from above are stored in a column-vector for easier manipulation (zbar).

The program calculates the weights that will yield a portfolio with minimum risk.

The risk of the portfolio is related to the standard deviation and it can be calculated as follows

So we now have to minimize:

Portfolio variance σ = √X'TΣX subject to: ∑iXi =1, where X is a vector of weights in portfolio
 
The expected return of the portfolio is given by w = X'(zbar)

The weights Xi are fractions invested in each security and thus they must add up to 1.

The minimization problem can be solved using Lagrange multipliers to minimize the variance.

### Calculate Efficient Frontier
This calculation uses the method outlines in Prof. Eric Zivot's Portfolio Theory Class notes here 

https://faculty.washington.edu/ezivot/econ424/portfolioTheoryMatrix.pdf

