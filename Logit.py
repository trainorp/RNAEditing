import numpy as np
import pandas as pd
import statsmodels.api as sm
# Fix missing chisq test
from scipy import stats
stats.chisqprob = lambda chisq, df: stats.chi2.sf(chisq, df)
from patsy import dmatrices

# Import data 
alt = pd.read_csv('/home/patrick/gdrive/RNAEditing/alt_matrix.chr1.tsv',delimiter='\t')
alt.shape
alt.dtypes
ref = pd.read_csv('/home/patrick/gdrive/RNAEditing/ref_matrix.chr1.tsv',delimiter='\t')
tot = pd.read_csv('/home/patrick/gdrive/RNAEditing/total_matrix.chr1.tsv',delimiter='\t')
ecpms = pd.read_csv('/home/patrick/gdrive/RNAEditing/ecpms.tsv',delimiter='\t')

# Group indicators:
g = range(2,23)
g1 = range(2,16)
g2 = range(16,23)
gInd = np.concatenate((np.repeat(1,len(g1)), np.repeat(0,len(g2))))

# One site example:
def siteLogisticTests(alt,ref,tot):
    resDf = pd.DataFrame({'chr': alt['chr'].values, 'pos': alt['pos'].values,
                      'diff': 0.0, 'p': 1.0, 'NT': ''})
    for i in range(alt.shape[0]):
        siteCntDf = pd.DataFrame({'n': tot.iloc[i,g],'Alt': alt.iloc[i,g],
                        'Ref': ref.iloc[i,g],'g': gInd, 'ecpms': ecpms['EPCM'].values})
    
        # Expanded site count DataFrame:
        siteCntExDf = pd.DataFrame(columns=['Alt','g'],dtype=int)
        for j in range(siteCntDf.shape[0]):
            # Binary indicators for alt vs ref:
            temp1 = np.concatenate((np.repeat(1,siteCntDf['Alt'].values[j]),
                            np.repeat(0,siteCntDf['Ref'].values[j])))
            # Repeat group indicator and EPCM
            temp2 = np.repeat(siteCntDf['g'].values[j],len(temp1))
            temp2b = np.repeat(siteCntDf['ecpms'].values[j],len(temp1))
            
            # Success/Fail & Group:
            temp3 = pd.DataFrame({'Alt': temp1, 'g': temp2, 'ecpms': temp2b})
            siteCntExDf = siteCntExDf.append(temp3)
        
        # All are in one group only: 
        if sum(siteCntExDf['g']) == siteCntExDf.shape[0]:
            resDf['NT'].values[i]='All in Ref Group'
            
        elif sum(siteCntExDf['g']) == 0:
            resDf['NT'].values[i]='All in Comparison Group'
            
        elif sum(siteCntExDf['Alt']) == siteCntExDf.shape[0]:
            resDf['NT'].values[i]='All Alt'
        
        # Else do chi-square test:
        else:
            # Logistic regression model:
            y, X = dmatrices('Alt~ecpms+g',data=siteCntExDf,return_type='dataframe')

            logit = sm.Logit(np.array(siteCntExDf['Alt'].values,dtype=int),X)
            
            # Fit model:
            try:
                logitRes = logit.fit() # Fit model
            except np.linalg.LinAlgError as e:
                resDf['NT'].values[i]=str(e)
            
            # Output:
            resDf['diff'].values[i] = logitRes.params['g[T.1]']
            resDf['p'].values[i] = logitRes.pvalues['g[T.1]']
            
    return resDf

# One site example:
def siteLogisticTests2(alt,ref,tot):
    resDf = pd.DataFrame({'chr': alt['chr'].values, 'pos': alt['pos'].values,
                      'diff': 0.0, 'p': 1.0, 'NT': ''})
    for i in range(alt.shape[0]):
        siteCntDf = pd.DataFrame({'n': tot.iloc[i,g],'Alt': alt.iloc[i,g],
                        'Ref': ref.iloc[i,g],'g': gInd, 'ecpms': ecpms['EPCM'].values,
                        'Prop': alt.iloc[i,g]/tot.iloc[i,g]})
        siteCntDf = siteCntDf[np.invert(np.isnan(siteCntDf['Prop'].values))]
        
        # Fix problems here:
        if len(set(siteCntDf['g'])) == 1:
            resDf['NT'].values[i]='All Alt in Same Group'
        
        elif len(set(siteCntDf['Alt'])) == 1:
            resDf['NT'].values[i]='No var in Alt/Ref'
            
        # Else do chi-square test:
        else:
            y, X = dmatrices('Prop~ecpms+g',data=siteCntDf,return_type='dataframe')
            logit = sm.GLM(y, X, family=sm.families.Binomial())
            
            try:
                logitRes = logit.fit()
            except sm.tools.sm_exceptions.PerfectSeparationError as e:
                resDf['NT'].values[i]=str(e)
                
            # Output:
            resDf['diff'].values[i] = logitRes.params['g']
            resDf['p'].values[i] = logitRes.pvalues['g']
            
    return resDf

siteLogisticTests(alt,ref,tot)