import numpy as np
import statsmodels.api as sm

# Fix missing chisq test
from scipy import stats
stats.chisqprob = lambda chisq, df: stats.chi2.sf(chisq, df)

# Fake data:
yUn = np.array([0,1,2,2,1,2,1,10,11,20,13,24,12,23,13,19])
w = np.array([26,27,28,29,30,31,22,23,24,25,26,27,28,29,30,31])
yEd = w - yUn

X = np.array([1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0])

# Expand data:
v1 = np.concatenate([np.repeat(1,yUn[x]) for x in range(0,len(yUn))])
v2 = np.concatenate([np.repeat(X[x],yUn[x]) for x in range(0,len(yUn))])
v = np.vstack((v1,v2)).T

v1 = np.concatenate([np.repeat(0,yEd[x]) for x in range(0,len(yEd))])
v2 = np.concatenate([np.repeat(X[x],yEd[x]) for x in range(0,len(yEd))])
v = np.vstack((v,np.vstack((v1,v2)).T))

# Statsmodels version
logit = sm.Logit(v[:,1],v[:,0]) # Model formula
logitRes = logit.fit() # Fit model
print(logitRes.summary()) # Print model 

logitRes.bse # Standard errors
logitRes.__dir__() # List all attributes
logitRes.pvalues # p-values