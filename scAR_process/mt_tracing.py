# get mt snp from bam file

import pypiper
outfolder = "hello_pypiper_results" # Choose a folder for your results

# Create a PipelineManager, the workhorse of pypiper
pm = pypiper.PipelineManager(name="hello_pypiper", outfolder=outfolder)

# Timestamps to delineate pipeline sections are easy:
pm.timestamp("Hello!")

# Now build a command and pass it to pm.run()
target_file = "hello_pypiper_results/output.txt"
command = "echo 'Hello, Pypiper!' > " + target_file
pm.run(command, target_file)

pm.stop_pipeline()



#
outfolder
samtools


import pypiper
import os

pm = pypiper.PipelineManager(name="mt_tracing", outfolder=outfolder)

# step 5 
# count number of reads on mt q30 PE
map_folder = os.path.join(outfolder, "Mapping")
ngstk.make_dir(map_folder)

mt_bam = os.path.join(map_folder, sample_name + ".mtQ30PE.bam")
cmd = tools.samtools + " view -b -q  30 " + " "
cmd += " -f 0x2 " + rmdup_bam + " "
cmd += " -o " + mt_bam + " "
cmd += " chrMT "
pm.run(cmd,mt_bam)

# count mt Q30PE
ar = ngstk.count_mapped_reads(mt_bam, args.paired_end)
pm.report_result("mtQ30PE_reads", ar)
#        try:
                # wrapped in try block in case Trimmed_reads is not reported in this pipeline.
tr = float(pm.get_stat("Aligned_reads"))
if tr!= 0: 
	pm.report_result("mtQ30PE_rate", round(float(ar) / float(tr), 2))
else:
	
	pm.report_result("mtQ30PE_rate", "0")
#        except:

#                pass

#  make re-align target regions  --step 5
#  Need to make genome wide re-align 
target_interval = os.path.join(map_folder, args.sample_name + ".target.intervals")
#Indel_vcf_s = os.path.join(args.output_parent, "mergeAll.mt.vcf")
#Indel_vcf_a = res.MTIndel           
cmd = "java -Xmx8g -jar " + tools.GATK + " "
cmd += " -R " + res.ref_genome_fasta + " " 
cmd += " -T RealignerTargetCreator "
#cmd += " -known " + Indel_vcf_s + " " 
cmd += " -known " + res.MTIndel + " "
cmd += " -known " + res.GIndel + " "
cmd += " -nt " + pm.cores + " " 
cmd += " -I " + rmdup_bam + " "
cmd += " -o " + target_interval
pm.run(cmd,target_interval)

# ReAlignment  --step 5 
realigned_bam = os.path.join(map_folder, args.sample_name + ".rmdup.realigned.bam")
cmd = "java -Xmx10g -jar " + tools.GATK + " "
cmd += "-R " + res.ref_genome_fasta + " " 
cmd += "-T IndelRealigner "
#cmd += " -known " + Indel_vcf_s + " "
#cmd += " -known " + Indel_vcf_a + " "
cmd += " -known " + res.GIndel + " "
cmd += " -known " + res.MTIndel + " "
#cmd += " -nt " + pm.cores + " "
cmd += " -I " + rmdup_bam + " "
cmd += " -targetIntervals " + target_interval
cmd += " -o " + realigned_bam
pm.run(cmd, realigned_bam)

count_steps += 1
if count_steps == args.stopN:
        print("pipeline stop at step")
        print(args.stopN)
        pm.stop


# mpileup 
#samtools mpileup  -l ../../chrM -f /home/jinxu/DB/hg19_g1k_v37/human_g1k_v37.fasta -q 20 -Q 20 13-TAAGGCGA-GCGATCTA.rmdup.realigned.recali.bam
mpileup = os.path.join(map_folder, args.sample_name + ".mpileup")
#chrMT =  res.chrM_len
cmd = tools.samtools + " mpileup " 
cmd += " -l " +  res.chrM_len + " "
cmd += " -f " + res.ref_genome_fasta + " "
cmd += " -q " + str(param.mpileup.q) + " "
cmd += " -Q " + str(param.mpileup.Q) + " "
cmd += param.mpileup.overlap + " "
cmd += realigned_bam
cmd += " > " + mpileup
pm.run(cmd,mpileup) 

# make count table for each allele 
count = os.path.join(map_folder, args.sample_name + ".counts")
cmd = tools.pileup_rj + " "
cmd += mpileup + " > "
cmd += count
pm.run(cmd,count)
#
# caculate average depth 
# awk 'BEGIN{sum=0;num=0}''{sum+=$4;num++}''END{print sum/num}'  count  
cmd = " awk 'BEGIN{sum=0;}''{sum+=$4;}''END{print sum}' " + count	
mt_total_depth = pm.checkprint(cmd)
cmd = " awk 'BEGIN{num=0;}''{num++;}''END{print num}' " + count	
mt_total_site =pm.checkprint(cmd)
#print(int(mt_total_site))
if int(mt_total_site)==0:
	pm.report_result("mt_aver_depth", "0")
	print(mt_total_site)
else:
	pm.report_result("mt_aver_depth", round(float(mt_total_depth) / float(mt_total_site), 2))
# call somatic mutation for each one  --step7
# java -jar /home/jinxu/software/VarScan.v2.3.7.jar   pileup2snp   13-TAAGGCGA-GCGATCTA.mpileup  --min-var-freq 0.001 --min-reads2 4
if int(mt_total_site) != 0:
	snv = os.path.join(map_folder, args.sample_name + ".snv")
	cmd = "java -Xmx4g  -jar " + tools.VarScan + " "
	cmd += " pileup2snp " +  mpileup + " "
	cmd += " --min-var-freq " + str(param.varscan.min_freq) + " " 
	cmd += " --min-reads2 " + str(param.varscan.min_r2) + " "
	cmd += " > " + snv
	pm.run(cmd,snv)
else:

	snv = os.path.join(map_folder, args.sample_name + ".snv")
# filter somatic mutation  --step6
# summarINPUT=test.samy using outside script.   --step7
pm.stop_pipeline()