# Exploring Relationships among International Stock Markets

In this project, we explore and analyze the performance of four international stock market indices based on their closing prices and net returns across time. The four selected international market indices are US S&P 500 (GSPC), Euro STOXX (STOXX), Hong Kong Hang Seng index (HSI) and Japan Nikkei 225 (N225). The sample time frame is between Jan 1, 2006 and Apr 30, 2021. The daily closing prices of the four selected international stock market indices, GSPC, STOXX, HSI,
N225, are extracted from [Yahoo! Finance](https://finance.yahoo.com/) by R package [“quantmod”](https://cran.r-project.org/web/packages/quantmod/quantmod.pdf).

The project is mainly divided into three parts. Here the R code accounts for the distribution analysis part. In this part, we explored the appropriate multivariate distribution and copula to four market indices’ net returns. We also performed goodness of fit.

The distribution analysis is divided into three parts:

### Part 1 
Univariate Distribution Analysis
- Perform Shapiro-Wilks normality test and conclude that all of the four returns do not follow the normal distribution
- Apply the Kolmogorov–Smirnov test, with the null hypothesis of sampling from the student t distribution and conclude that these four indices' returns may follow the univariate t distribution
- We could obtain the fitted MLE “std” t distribution parameters by applying R function “fitdist()”

### Part 2
Multivariate Distribution Analysis
- Fit the multivariate t distribution by the profile likelihood method to obtain the fitted "MLE" parameters in multi-dimension, i.e. mean vector μ, scale matrix Λ and tail index ν
- We repeatedly call the R function “cov.trob()” function to find the optimal value of the tail index and determine the MLE of the other two candidates

### Part 3
T-copula and Meta-t Distribution
- Determine Kendall's tau correlation coefficients and transform into the correlation matrix of our net return matrix
- Transform the net returns into a univariate uniform and apply the R functions “tCopula()” and “fitCopula()” to fit the t-copula to our four indices' return matrix
- Define the meta-t distribution to our net return matrix given the “std” marginal distribution by applying R functions “mvdc()” and “fitMvdc()
- According to the fitted meta-t distribution parameters, GSPC and STOXX’s net returns are highly correlated. Besides, HSI and N225’s net returns also have high
correlation, which may indicate the possibility of “co-movement” in their underlying closing prices.

Conclusion: In distribution analysis, we found that the US and the European markets show strong correlation and interdependence, as well as the Hong Kong and the Japanese markets when analyzing their returns performance. This project shows the possibility of applying cointegration technique when constructing portfolios with these four indices to induce stationary time series process, but more practical researches are needed for further analysis.
