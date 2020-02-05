'''
inputdir='/Share2/home/zhangqf5/yanqiu/scAR/20190925/'
outf='log/20190925_ATAC_raw.list'
with open(outf,'w') as output:
	for i in range(9):
		line="""A{0}\tTL000A{0}.\t{1}/TL_A{0}_1.fq.gz\t{1}/TL_A{0}_2.fq.gz""".format(i+1, inputdir)
		output.write(line+'\n')

inputdir='/Share2/home/zhangqf5/yanqiu/scAR/20190925/'
outf='log/ZJS0925_A2_raw.list'
with open(outf,'w') as output:
	for i in range(9):
		line="""A{0}\tZJS00A{0}.\t{1}/ZJS2_A{0}_1.fq.gz\t{1}/ZJS2_A{0}_2.fq.gz""".format(i+1, inputdir)
		output.write(line+'\n')
	for i in range(9,24):
		line="""A{0}\tZJS0A{0}.\t{1}/ZJS2_A{0}_1.fq.gz\t{1}/ZJS2_A{0}_2.fq.gz""".format(i+1, inputdir)
		output.write(line+'\n')


outf='log/ZJS0925_A1_raw.list'
with open(outf,'w') as output:
	for i in range(9): 
		line="""A{0}\tZJS00A{0}.\t{1}/ZJS1_A{0}_1.fq.gz\t{1}/ZJS1_A{0}_2.fq.gz""".format(i+1, inputdir)
		output.write(line+'\n')

inputdir='/Share2/home/zhangqf5/yanqiu/scAR/ZJSMix0928/'
outf='log/ZJSMix0928_coA_raw.list'
with open(outf,'w') as output:
	for i in range(9):
		line="""coA{0}\tCO000A{0}.\t{1}/coA{0}_1.fq.gz\t{1}/coA{0}_2.fq.gz""".format(i+1, inputdir)
		output.write(line+'\n')
	line="""coA{0}\tCO00A{0}.\t{1}/coA{0}_1.fq.gz\t{1}/coA{0}_2.fq.gz""".format(10, inputdir)
	output.write(line)


inputdir='/Share2/home/zhangqf5/yanqiu/scAR/20191025'
outf='log/ZJS191025_NorR_raw.list'
with open(outf,'w') as output:
	for i in range(9):
		line="""Nor{0}\tZJSNor{0}.\t{1}/Nor_R{0}_1.fq.gz\t{1}/Nor_R{0}_2.fq.gz""".format(i+1, inputdir)
		output.write(line+'\n')
	line="""Nor{0}\tZJSNo{0}.\t{1}/Nor_R{0}_1.fq.gz\t{1}/Nor_R{0}_2.fq.gz""".format(10, inputdir)
	output.write(line)


inputdir='/Share2/home/zhangqf5/yanqiu/scAR/20191030'
outf='log/ZJS191030_A_raw.list'
with open(outf,'w') as output:
	for i in range(3):
		line="""N{0}\tZJS00N{0}.\t{1}/ZJS_N{0}A_1.fq.gz\t{1}/ZJS_N{0}A_2.fq.gz""".format(i+1, inputdir)
		output.write(line+'\n')

inputdir='/Share2/home/zhangqf5/yanqiu/scAR/20191116'
outf='log/ZJS191116_coA_raw.list'
with open(outf,'w') as output:
	for i in range(8):
		line="""coA{0}\tZJS0co{0}.\t{1}/ZJScoA{0}_1.fq.gz\t{1}/ZJScoA{0}_2.fq.gz""".format(i+1, inputdir)
		output.write(line+'\n')

inputdir='/Share2/home/zhangqf5/yanqiu/scAR/20191116'
outf='log/ZJS191116_coR_raw.list'
with open(outf,'w') as output:
	for i in range(8):
		line="""coR{0}\tZJS0co{0}.\t{1}/ZJScoR{0}_1.fq.gz\t{1}/ZJScoR{0}_2.fq.gz""".format(i+1, inputdir)
		output.write(line+'\n')

outf='log/ZJS191116_oR_raw.list'
with open(outf,'w') as output:
	i=8
	line="""oR{0}\tZJS00o{0}.\t{1}/ZJSoR{0}_1.fq.gz\t{1}/ZJSoR{0}_2.fq.gz""".format(i+1, inputdir)
	output.write(line+'\n') 
	for i in range(9,16):
		line="""oR{0}\tZJS0o{0}.\t{1}/ZJSoR{0}_1.fq.gz\t{1}/ZJSoR{0}_2.fq.gz""".format(i+1, inputdir)
		output.write(line+'\n')
outf='log/ZJS191116_oA_raw.list'
with open(outf,'w') as output:
	for i in range(16,24):
		line="""oA{0}\tZJS0o{0}.\t{1}/ZJSoA{0}_1.fq.gz\t{1}/ZJSoA{0}_2.fq.gz""".format(i+1, inputdir)
		output.write(line+'\n')

outf='log/ZJS200110_A_raw.list'
inputdir='/Share2/home/zhangqf5/yanqiu/scAR/20200110'
with open(outf,'w') as output:
	for i in [1,8]:
		line="""Z{0}\tZJS000{0}.\t{1}/ZJS_{0}A_1.fq.gz\t{1}/ZJS_{0}A_2.fq.gz""".format(i, inputdir)
		output.write(line+'\n')

outf='log/ZJS200110_R_raw.list'
with open(outf,'w') as output:
    for i in [1,8]:
        line="""Z{0}\tZJS000{0}.\t{1}/ZJS_{0}R_1.fq.gz\t{1}/ZJS_{0}R_2.fq.gz""".format(i, inputdir)
        output.write(line+'\n') 
'''


inputdir='/Share2/home/zhangqf5/yanqiu/scAR/20200110'
outf='log/WYF200110_12A_raw.list'
with open(outf,'w') as output:
	for i in [12,34]:
		line="""W{0}\tWYF00{0}.\t{1}/WYF_{0}A_1.fq.gz\t{1}/WYF_{0}A_2.fq.gz""".format(i, inputdir)
		output.write(line+'\n')
outf='log/WYF200110_12R_raw.list'
with open(outf,'w') as output: 
	for i in [12,34]:
		line="""W{0}\tWYF00{0}.\t{1}/WYF_{0}R_1.fq.gz\t{1}/WYF_{0}R_2.fq.gz""".format(i, inputdir)
		output.write(line+'\n')

