from plot_cor import plot_cor
import pandas as pd

'''
# generate gene annotation
import pickle
def get_gene_anno(gene_info):
    gene_id_to_name=gene_info['gene_id_to_name']
    gene_name_to_id={}
    for gene_id in gene_id_to_name.keys():
        gene_name=gene_id_to_name[gene_id]
        gene_name_to_id.setdefault(gene_name,[])    
        gene_name_to_id[gene_name].append(gene_id)
        gene_name_to_id.setdefault(gene_id.split('.')[0],[])
        gene_name_to_id[gene_id.split('.')[0]].append(gene_id)
        gene_name_to_id.setdefault(gene_id,[])
        gene_name_to_id[gene_id].append(gene_id)
    gene_anno={'id_to_name':gene_info['gene_id_to_name'], 'name_to_id':gene_name_to_id}    
    return gene_anno

gene_info_hg=pickle.load(open('/Share2/home/zhangqf5/yanqiu/scAR_old/split_seq/split-seq-pipeline/INDEX_hg/gene_info.pkl','rb'))
gene_info_mm=pickle.load(open('/Share2/home/zhangqf5/yanqiu/scAR_old/split_seq/split-seq-pipeline/INDEX_mm/gene_info.pkl','rb'))
gene_anno_hg=get_gene_anno(gene_info_hg)
gene_anno_mm=get_gene_anno(gene_info_mm)
pickle.dump(gene_anno_hg,open('/Share2/home/zhangqf5/yanqiu/library/gene_anno_hg.pkl','wb') )
pickle.dump(gene_anno_mm, open('/Share2/home/zhangqf5/yanqiu/library/gene_anno_mm.pkl','wb'))
'''

def annotate_gene(genes, gene_anno):
    filtered=[]
    gene_ids=[] 
    for gene in genes:
        try:
            gene_id=gene_anno['name_to_id'][gene]
            if len(gene_id)==1 and (gene_id[0] not in gene_ids):
                filtered.append(gene)
                gene_ids.append(gene_id[0])
        except:
            pass
    return filtered, gene_ids

def is_gene_id(name):
    if name.startswith('ENS') and len(name.split('.'))==2:
        return True

def plot_cor_df(fpkm1, fpkm2, sample1,sample2,gene_anno,normalize=1000000, how='inner',with_density=False, save=None):
    fpkm1=pd.read_csv(fpkm1, sep='\t', header=0, index_col=0)
    fpkm2=pd.read_csv(fpkm2, sep='\t', header=0, index_col=0)
    fpkm1 = fpkm1.loc[~fpkm1.index.duplicated(keep='first')]
    fpkm2 = fpkm2.loc[~fpkm2.index.duplicated(keep='first')]
    if not is_gene_id(fpkm1.index.values[0]):
        filtered, gene_ids=annotate_gene(fpkm1.index.values, gene_anno)
        fpkm1=fpkm1.loc[filtered,:]
        fpkm1.index=gene_ids
    if not is_gene_id(fpkm2.index.values[0]):
        filtered, gene_ids=annotate_gene(fpkm2.index.values, gene_anno)
        fpkm2=fpkm2.loc[filtered,:]
        fpkm2.index=gene_ids
    fpkm1.index=[x.split('.')[0] for x in fpkm1.index]
    fpkm2.index=[x.split('.')[0] for x in fpkm2.index]
    merged=pd.merge(fpkm1, fpkm2, left_index=True, right_index=True, how=how)
    merged.columns=[sample1,sample2 ]
    cor=plot_cor(merged[sample1], merged[sample2],sample1,sample2, normalize=normalize,with_density=with_density, save=save)
    
if __name__=='__main__':
    import sys
    fpkm1=sys.argv[1]
    fpkm2=sys.argv[2]
    sample1=sys.argv[3]
    sample2=sys.argv[4]
    genome=sys.argv[5] # for annotating gene name, mm or hg
    save=sys.argv[6]
    
    import pickle
    with open('/Share2/home/zhangqf5/yanqiu/scAR_old/split_seq/split-seq-pipeline/INDEX_%s/gene_info.pkl'%genome, 'rb') as f:
        gene_info = pickle.load(f)

    with open('/Share2/home/zhangqf5/yanqiu/library/gene_anno_%s.pkl'%genome, 'rb') as f:
        gene_anno = pickle.load(f)
        
    plot_cor_df(fpkm1, fpkm2, sample1,sample2, gene_anno, normalize=1000000, how='inner',with_density=False, save=save)    
    