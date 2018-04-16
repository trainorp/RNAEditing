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
def siteTests(alt,ref,tot):
    resDf = pd.DataFrame({'chr': alt['chr'].values, 'pos': alt['pos'].values,
                      'diff': 0.0, 'p': 1.0, 'NT': ''})
    for i in range(alt.shape[0]):
        siteCntDf = pd.DataFrame({'n': tot.iloc[i,g],'Alt': alt.iloc[i,g],
                                  'Ref': ref.iloc[i,g],'g': gInd})
    
        # Expanded site count DataFrame
        siteCntExDf = pd.DataFrame(columns=['Alt','g'],dtype=int)
        
        for j in range(siteCntDf.shape[0]):
            temp1 = np.concatenate((np.repeat(1,siteCntDf['Alt'].values[j]),
                            np.repeat(0,siteCntDf['Ref'].values[j])))
            temp2 = np.repeat(siteCntDf['g'].values[j],len(temp1))
            temp3 = pd.DataFrame(np.stack((temp1,temp2)).T,columns=['Alt','g'])
            siteCntExDf = siteCntExDf.append(temp3)
            
            if sum(siteCntExDf['g']) == siteCntExDf.shape[0]:
                resDf['NT'].values[i]='a'
                
            elif sum(siteCntExDf['g']) == 0:
                resDf['NT'].values[i]='b'
                
            else:
                # Logistic regression model:
                logit = sm.Logit(np.array(siteCntExDf.iloc[:,1],dtype=int),
                                 np.array(siteCntExDf.iloc[:,0],dtype=int))
                logitRes = logit.fit() # Fit model
                
                # Output:
                resDf['diff'].values[i] = logitRes.params[0]
                resDf['p'].values[i] = logitRes.pvalues
        