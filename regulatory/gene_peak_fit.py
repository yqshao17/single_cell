import sys
import numpy as np
import pandas as pd
import scanpy as sc
import pickle

    
def get_peak(peak):
    if ':' in peak: # chr:start-end
        ch=peak.split(':')[0]
        s,e=peak.split(':')[1].split('-')
    elif '_' in peak: # chr_start_end
        ch,s,e=peak.split('_')
    return ch,int(s),int(e)

def intersect_features(features, genome_region):
    chosen=[]
    for feature in features:
        ch,s,e=get_peak(feature)
        if ch==genome_region[0]:
            if genome_region[1]<e<genome_region[2] or genome_region[1]<s<genome_region[2]:
                chosen.append(feature)
    return chosen

# linear regression
import statsmodels.api as sm
import scipy.stats

def get_gene_peak_ols(gene, around_peaks, rn, an):
    atac_df={}
    for peak in around_peaks:
        atac_df[peak]=an.obs_vector(peak)
    X = pd.DataFrame(atac_df)
    y = rn.obs_vector(gene)
    est = sm.OLS(y, X).fit()
    resid=est.resid
    basic=est.summary2().tables[0]
    coefs=est.summary2().tables[1]
    R2=float(basic.iloc[6,1])
    R2_adj=float(basic.iloc[0,3])
    coefs=coefs.astype('float')        
    filtered_peaks=coefs[coefs['P>|t|']<0.05].sort_values('Coef.', ascending=False)
    return R2, R2_adj, resid, filtered_peaks

def get_gene_peak_ols(gene, around_peaks, rn, an, weights):
    atac_df={}
    for peak in around_peaks:
        atac_df[peak]=an.obs_vector(peak)
    X = pd.DataFrame(atac_df)
    y = rn.obs_vector(gene)
    est = sm.WLS(y, X, weights=weights).fit()
    predictions = est.predict(X)
    resid=est.resid
    basic=est.summary2().tables[0]
    coefs=est.summary2().tables[1]
    R2=float(basic.iloc[6,1])
    R2_adj=float(basic.iloc[0,3])
    coefs=coefs.astype('float')
    peaks=coefs.sort_values('Coef.', ascending=False)
    filtered_peaks=coefs[coefs['P>|t|']<0.05].sort_values('Coef.', ascending=False)
    
    return R2, R2_adj, resid, peaks, filtered_peaks, predictions



def get_gene_peak_lasso(gene, around_peaks, rn, an):
    
    return
    
    

def get_gene_peak_cor(gene, around_peaks, rn, an, method='spearman'):
    def cal_cor(a, b):
        if method=='spearman':
            cor=scipy.stats.spearmanr(a,b) 
        elif method=='pearson':
            cor=scipy.stats.pearsonr(a,b) 
        return cor
    peaks=pd.DataFrame(index=['Coef.','pval'])
    filtered_peaks=pd.DataFrame(index=['Coef.','pval'])
    for peak in around_peaks:
        cor=cal_cor(an.obs_vector(peak),rn.obs_vector(gene))
        peaks[peak]=[cor[0],cor[1]]
        if cor[1]<0.05 or cor[0]>0.25:
            filtered_peaks[peak]=[cor[0],cor[1]]
    peaks=peaks.T
    filtered_peaks=filtered_peaks.T
    return peaks, filtered_peaks
        
        
def get_gene_peak_fit(gene_info, rn, an, outprefix, bin_size, method='regress'):
    # get gene-peak regulatory relationship
    bin_size=int(bin_size)
    gene_id_to_name=gene_info['gene_id_to_name']
    gene_name_to_id={}
    for gene_id in gene_id_to_name.keys():
        gene_name=gene_id_to_name[gene_id]
        gene_name_to_id.setdefault(gene_name,[])    
        gene_name_to_id[gene_name].append(gene_id)
        gene_name_to_id.setdefault(gene_id.split('.')[0],[])
        gene_name_to_id[gene_id.split('.')[0]].append(gene_id)
    gene_id_to_chrom=gene_info['gene_id_to_chrom']
    gene_starts=gene_info['gene_starts']
    gene_ends=gene_info['gene_ends']
    
        
    def get_gene_id(gene):
        if gene in gene_id_to_name.keys():
            gene_id = gene
        elif gene in gene_name_to_id.keys():
            gene_id = gene_name_to_id[gene][0]
            if len(gene_name_to_id[gene])>1:
                log.write(gene+' has ambiguous gene name. Default choose %s!\n'%gene_id)
        else:
            log.write(gene+' has no valid annotation!\n')
            raise NameError('Not valid gene name: %s!'%gene)
        return gene_id
    
    def get_around_peaks(gene_id, around_bin, peaks):       
        chrom=gene_id_to_chrom[gene_id].split('_')[-1]
        #chrom=gene_id_to_chrom[gene_id]
        start=gene_starts[gene_id]  
        end=gene_ends[gene_id]
        start_bin=start-around_bin
        end_bin=end+around_bin
        around_peaks=intersect_features(peaks, (chrom, max(0,start_bin), end_bin))
        return around_peaks
    
    # get peaks
    
    log=open(outprefix+'gene_around_peaks_%dk.log'%bin_size,'w')
    
    gene_peaks={}
    for gene in rn.var_names:
        gene_peaks[gene]=[]
        try:
            gene_id=get_gene_id(gene)
        except:
            pass
        peaks=get_around_peaks(gene_id, bin_size*1000, an.var_names)
        log.write('%s: %d around peaks\n'%(gene, len(peaks)))
        gene_peaks[gene]=peaks
            
    with open(outprefix+'gene_around_peaks_%dk.pkl'%bin_size, 'wb') as f:
        pickle.dump(gene_peaks, f)  
    log.close()
    
    
    #gene_peaks=pickle.load(open(outprefix+'gene_around_peaks_%dk.pkl'%bin_size,'rb'))
    
    if method=='regress':
        #gene_peaks = pickle.load(open(outprefix+'/gene_around_peaks_250k.pkl', 'rb'))
        gene_peak_fit={'peaks':{},'filtered_peaks':{},'R2':{},'R2_adj':{},'resid':{},'predictions':{}}
        log=open(outprefix+'gene_peak_regress_%dk.log'%bin_size,'w')
        weights=rn.obs['Weights'].values
        for gene in gene_peaks.keys():
            around_peaks=gene_peaks[gene]
            if around_peaks:
                R2, R2_adj, resid, peaks, filtered_peaks, predictions=get_gene_peak_ols(gene, around_peaks, rn, an, weights)
                gene_peak_fit['peaks'][gene]=peaks
                gene_peak_fit['filtered_peaks'][gene]=filtered_peaks
                gene_peak_fit['R2'][gene]=R2
                gene_peak_fit['R2_adj'][gene]=R2_adj
                gene_peak_fit['resid'][gene]=resid
                gene_peak_fit['predictions'][gene]=predictions
                log.write('%s: %d regulatory peaks\n'%(gene, filtered_peaks.shape[0]))
        log.close()  
    
    elif method=='spearman' or method=='pearson':
        # gene peak cor
        gene_peak_fit={'peaks':{},'filtered_peaks':{}}
        log=open(outprefix+'gene_peak_%s_%dk.log'%(method,bin_size),'w')
        for gene in gene_peaks.keys():
            around_peaks=gene_peaks[gene]
            if around_peaks:
                peaks, filtered_peaks=get_gene_peak_cor(gene,around_peaks,rn,an)
                gene_peak_fit['peaks'][gene]=peaks
                gene_peak_fit['filtered_peaks'][gene]=filtered_peaks
                log.write('%s: %d regulatory peaks\n'%(gene, filtered_peaks.shape[0]))
        log.close()
        
    with open(outprefix+'gene_peak_%s_%dk.pkl'%(method,bin_size), 'wb') as f:
        pickle.dump(gene_peak_fit, f)
        
    return gene_peak_fit

    
# get regulatory matrix

def get_obs_assign(gene, peak, coef,  rn, an): #1
        assignment=[]
        for exp, acc in zip(rn.obs_vector(gene), an.obs_vector(peak)):
            if exp>0 and acc>0:
                assignment.append(1)
            else:
                assignment.append(0)
        return assignment

'''
def get_obs_assign(gene, peak, coef, rn, an): #2
        assignment=[]
        for acc in an.obs_vector(peak):
            assignment.append(acc*coef)
        return assignment


def get_obs_assign(gene, peak, coef, rn, an): #3
        assignment=[]
        for exp, acc in zip(rn.obs_vector(gene), an.obs_vector(peak)):
            if exp>0:
                assignment.append(acc*coef/exp)
            else:
                assignment.append(0)
        return assignment    
'''

def get_regulatory_matrix(gene_peak_fit_coef, rn, an):   
    # gene_pair cell assignment
    # rn and an should have the same obs_names
    obs_assign={}
    for gene in gene_peak_fit_coef.keys():
        peakdf=gene_peak_fit_coef[gene]
        if not peakdf.empty:
            for peak in peakdf.index:
                coef=peakdf.loc[peak,'Coef.']
                obs_assign[gene+'_'+peak]=get_obs_assign(gene, peak, coef, rn, an)
    reg_df=pd.DataFrame(obs_assign)
    reg_df.index=rn.obs_names
    reg_df=reg_df.T
    return reg_df

def pseudo_regulatory_matrix(gene_peak_fit_coef, CellID):   #
    # gene_pair cell assignment
    # rn and an should have the same obs_names
    obs_assign={}
    for gene in gene_peak_fit_coef.keys():
        peakdf=gene_peak_fit_coef[gene]
        if not peakdf.empty:
            for peak in peakdf.index:
                coef=peakdf.loc[peak,'Coef.']
                obs_assign[gene+'_'+peak]=coef
    reg_df=pd.DataFrame()
    reg_df[CellID]=obs_assign.values()
    reg_df.index=obs_assign.keys()
    return reg_df
