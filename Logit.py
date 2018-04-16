import numpy as np
import pandas as pd
import statsmodels.api as sm

# Fix missing chisq test
from scipy import stats
stats.chisqprob = lambda chisq, df: stats.chi2.sf(chisq, df)

# Import data 
alt = pd.read_csv('/home/patrick/gdrive/RNAEditing/alt_matrix.chr1.tsv',delimiter='\t')
alt.shape
alt.dtypes
ref = pd.read_csv('/home/patrick/gdrive/RNAEditing/ref_matrix.chr1.tsv',delimiter='\t')
tot = pd.read_csv('/home/patrick/gdrive/RNAEditing/total_matrix.chr1.tsv',delimiter='\t')

# Group indicators:
g = range(2,23)
g1 = range(2,16)
g2 = range(16,23)
gInd = np.concatenate((np.repeat(1,len(g1)), np.repeat(0,len(g2))))

# One site example:
i = 0 
siteCntDf = pd.DataFrame({'n': tot.iloc[i,g],'Alt': alt.iloc[i,g],'Ref': ref.iloc[i,g],
                            'g': gInd})
siteCntExDf = pd.DataFrame(columns=['Alt','g'],dtype=int)
for j in range(siteCntDf.shape[0]):
    temp1 = np.concatenate((np.repeat(1,siteCntDf['Alt'].values[j]),
                    np.repeat(0,siteCntDf['Ref'].values[j])))
    temp2 = np.repeat(siteCntDf['g'].values[j],len(temp1))
    temp3 = pd.DataFrame(np.stack((temp1,temp2)).T,columns=['Alt','g'])
    siteCntExDf = siteCntExDf.append(temp3)
    
logit = sm.Logit(np.array(siteCntExDf.iloc[:,1],dtype=int),
                 np.array(siteCntExDf.iloc[:,0],dtype=int))
logitRes = logit.fit() # Fit model
print(logitRes.summary()) # Print model 

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