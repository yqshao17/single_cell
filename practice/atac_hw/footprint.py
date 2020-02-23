import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from optparse import OptionParser
from optparse import OptionGroup
import numpy as np
from utils import Mkdir
import os

def Bed2bedGraph(Inbed,ref_size):
    Outbedgraph=Inbed[:-4]+'.bedGraph'
    os.system("genomeCoverageBed -bg -split -i %s -g %s > %s" %(Inbed,ref_size,Outbedgraph))
    return

def PerbasebedGraph(InbedGraph):
    tmp=InbedGraph+'.tmp'
    tmpf=open(tmp,'w')
    lines=open(InbedGraph,'r')
    for line in lines:
        items=line.split('\t')
        for i in range(int(items[2])-int(items[1])):
            tmpf.write(items[0]+'\t'+str(int(items[1])+i)+'\t'+str(int(items[1])+i+1)+'\t'+str(items[3]))
    os.system('mv %s %s'%(tmp,InbedGraph))
    return   

def Drawfootprint(incsv,motifname,halfwidth):
    lines = open(incsv).readlines()
    length=len(lines[0].split('\t'))
    counts_matrix, genes = [], []
    for iline,line in enumerate(lines):
        if iline>0:
            words = line.split('\t')
            #counts_matrix.append(map(float, words[1:]))
            counts_matrix.append([float(el) for el in words[1:]])
            genes.append(words[0])
    counts_matrix = np.asarray(counts_matrix)
    n_ave, n_step, cut_off = 50, 10, 0.7
    n_limit = (len(counts_matrix) - n_ave) // n_step
    ave_matrix = [counts_matrix[i*n_step:i*n_step+n_ave, :].sum(axis=0)/float(n_ave) for i in range(0, n_limit)]
    ave_matrix = np.asarray(ave_matrix)
    log_ave_matrix = np.log(ave_matrix+1)
    max_logAve = log_ave_matrix.max()
    log_ave_matrix[np.where(log_ave_matrix > max_logAve*cut_off)] = max_logAve*cut_off
    nt_label = [' '] * (int(length)-1)
    nt_label[0], nt_label[int(len(nt_label)/2)],nt_label[-1] = '-'+str(halfwidth),motifname,str(halfwidth)
    log_ave_df = pd.DataFrame(log_ave_matrix, index=range(0, len(log_ave_matrix)), columns=nt_label)
    sns.set_context('poster', font_scale=1)
    fig=plt.figure(figsize=(6, 14))
    ax11 = plt.subplot2grid((60, 1), (14, 0), rowspan=46)
    sns.heatmap(log_ave_df, xticklabels=True, yticklabels=False, cbar=False, ax=ax11)
    plt.setp(ax11, ylabel=motifname+' motif-sites')
    ax12 = plt.subplot2grid((60, 1), (0, 0), rowspan=9)
    ax12.plot(range(-(halfwidth), halfwidth+1), log_ave_df.sum(axis=0), color='firebrick')
    ax12.set_xlim([-(halfwidth), (halfwidth)])
    ax12.set_xticks([-halfwidth, 0, halfwidth])
    ax12.set_xlabel('Distance to motif')
    ax12.set_ylabel('Insertion-site probability')
    ax12.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
    fig.savefig(incsv+'.png', bbox_inches='tight')
    plt.close(fig)
    
    
#######  Options  #######

usage="usage: %prog [options][inputs]"
parser=OptionParser(usage=usage, version="%prog 1.0")

parser.add_option("--ref_size",help="reference chrom_size")
parser.add_option("-r",help="reference name")
parser.add_option("--motif",help="set motif")
#parser.add_option("--fragfile",help="fragment file")
parser.add_option("--perbase",help="set perbase bed file")
parser.add_option("--peak",help="peak file")
parser.add_option("--window",default=250,help="set window width,default=500")
parser.add_option("--norm",default=5000000, help="set norm num,default=5000000")
parser.add_option("-o",type="string",default='output',help="Set output directory")
parser.add_option("-t",type='int',default=4,help="Set the max thread,default = 4")


(options, args) = parser.parse_args()

def Footprint():
    peaklist=options.peak
    ref=options.r
    inmotif=options.motif
    outfootprint=options.o

    #Mkdir(outfootprint)
    bg=options.perbase[:-4]+'.bedGraph'
    if not os.path.exists(bg):
        Bed2bedGraph(options.perbase,options.ref_size)
        PerbasebedGraph(bg)
    imname='/home/xionglei/yanqiu/lib/Homer_motif/'+inmotif+'.motif'
    #mlen=len(open(imname).readlines())-1
    window=int(options.window)
    omname=outfootprint+'/'+inmotif
    peakname=peaklist.split('/')[-1]
    opeakname=outfootprint+'/'+peakname
    outcount=opeakname+'.'+inmotif+'_'+bg.split('/')[-1]+'.count'
    outvalue=opeakname+'.'+inmotif+'_'+bg.split('/')[-1]+'.value'
    outtxt=opeakname+'.'+inmotif+'_'+bg.split('/')[-1]+'.footprint'
    
    if (not os.path.exists(imname+'.'+ref+'.bed') or os.path.getsize(imname+'.'+ref+'.bed')==0):
        CMD = "scanMotifGenomeWide.pl "+imname+' '+ref+' -bed > '+imname+'.'+ref+'.bed'
        print('CMDLINE:'+CMD)
        os.system(CMD)
    if (not os.path.exists(opeakname+'.'+inmotif+'.in') or os.path.getsize(opeakname+'.'+inmotif+'.in')==0):
        os.system("bedtools intersect -a "+imname+'.'+ref+'.bed'+' -b '+peaklist+' -wo > '+opeakname+'.'+inmotif+'.in')

    if (not os.path.exists(opeakname+'.'+inmotif+'.position') or os.path.getsize(opeakname+'.'+inmotif+'.position')==0):
        inlist=pd.read_table(opeakname+'.'+inmotif+'.in',header=None)
        center=np.array((inlist[1]+inlist[2])/2).astype('int')
        inlist[1]=center-window
        inlist[2]=center+window+1
        inlist.loc[:,0:2].to_csv(opeakname+'.'+inmotif+'.position',header=False,sep='\t',index=False)
    if (not os.path.exists(outcount) or os.path.getsize(outcount)==0):    
        os.system("bedtools intersect -a "+bg+" -b "+opeakname+'.'+inmotif+'.position -wa -wb > ' + outcount)
    if (not os.path.exists(outvalue) or os.path.getsize(outvalue)==0):
        a=pd.read_table(outcount,header=None,index_col=None)
        a['center']=((a[5]+a[6])/2).astype(int)
        a['dis']=a[1]-a['center']
        coj=lambda x: ''.join([x[0],':',str(x['center'])])
        a['pos']=a.apply(coj,axis=1)
        a['value']=a[3]
        a[['pos','dis','value',]].to_csv(outvalue,header=True,index=False,sep='\t')
        b=a[['pos','dis','value',]]
    else:
        b=pd.read_table(outvalue,header=0,index_col=None)

    df=pd.DataFrame(np.zeros((len(set(b['pos'])), int(window)*2+1)),index=set(b['pos']) )
    df.columns=range(-int(window),int(window)+1)
    for i in range(b.shape[0]):
        df.loc[b.iloc[i,0],b.iloc[i,1]]=b.iloc[i,2]

    mysum=int(os.popen('wc -l %s'%(options.perbase)).read().split()[0])

    df['sum']=df.sum(axis=1)
    df=df.sort_values(by=["sum"],ascending=False).drop(['sum'],axis=1)
    if int(options.norm)!=0:
        df=df*int(options.norm)/mysum

    df.to_csv(outtxt,sep='\t',index=True, header=True)
    df.sum().to_csv(outtxt+'.sum',sep='\t',index=True)
    
    #Drawfootprint(outtxt,inmotif,window)

Footprint()
    
'''
ref='hg19'
ref_size='/home/xionglei/yanqiu/lib/hg19.chrom.sizes'
outfootprint='result/Hematopoiesis-All/footprint'

perbase='result/Hematopoiesis-All/footprint/Cluster3_fragment.bed'
bg=perbase[:-4]+'.bedGraph'
window=500
inmotif='gata2' # motif from homer database http://homer.ucsd.edu/homer/custom.motifs
omname=outfootprint+'/'+inmotif
peaklist='data/Hematopoiesis-All/Hematopoiesis-All_peak.bed'
peakname=peaklist.split('/')[-1]
opeakname=outfootprint+'/'+peakname
imname='/home/xionglei/yanqiu/lib/Homer_motif/'+inmotif+'.motif'
outtxt=opeakname+'.'+inmotif+'_'+bg.split('/')[-1]+'.footprint'
'''

#### usage ###
'''
export PATH=~/yanqiu/package/homer/bin:$PATH
python footprint.py \
    -r hg19 --ref_size /home/xionglei/yanqiu/lib/hg19.chrom.sizes \
    -o result/Hematopoiesis-All/footprint \
    --peak data/Hematopoiesis-All/Hematopoiesis-All_peak.bed \
    --motif ebf \
    --perbase result/Hematopoiesis-All/footprint/Cluster2_fragment.bed 
'''    