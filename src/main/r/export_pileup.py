import sys
import pysam

def export_sample(sample, ref):
	infile = "bam/"+sample+"."+ref+".bam"
	outfile = "variants/"+sample+"."+ref+".txt" 
	out = open(outfile, "w")
	out.write('position\tread\tnt\n')
	samfile = pysam.Samfile(infile, "rb" )
	for pileupcolumn in samfile.pileup( ref):
	    position =  pileupcolumn.pos
	    #print 'coverage at base %s = %s' % (pileupcolumn.pos , pileupcolumn.n)
	    for pileupread in pileupcolumn.pileups:
		qname = pileupread.alignment.qname
		#arr = qname.split(":")
		#read = arr[2] + ":" + arr[3]
		read = qname
	    	nt = pileupread.alignment.seq[pileupread.qpos]
		#print '\tbase in read %s = %s' % (pileupread.alignment.qname, pileupread.alignment.seq[pileupread.qpos])
		line ='%s\t%s\t%s\n' % (position, read, nt)
		#print line
		out.write(line)

	out.close()
	samfile.close()

def export_samples2():
	#for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]:
	for i in range(20):
		num = i+1
		if num < 10: num='0%s' % num
		sample = 'sample%s' % num
		print sample

def export_samples(samples, refs):
	for num in samples:
		if num < 10: num='0%s' % num
		sample = 'sample%s' % num
		for ref in refs:
			print '%s.%s' % (sample, ref)
			export_sample(sample,ref)

# python export_pileup.py sample18 HHaa36
sample = sys.argv[1] 
ref = sys.argv[2]

export_sample(sample,ref)

#hh = ['HHaa36','HHaa156']
#nk = ['NKaa36']
#kt9 = ['KT9aa36','KT9aa156']

#export_samples([1,2,3], hh)
#export_samples([4,5,6,7], nk)
#export_samples([8,9,10,11,12,13,14,15,16,17], kt9)
#export_samples([18,19,20], hh)

#hhvar = ['HHaa36var','HHaa156var']
#nkvar = ['NKaa36var']
#kt9var = ['KT9aa36var','KT9aa156var']

#export_samples([1,2,3], hhvar)
#export_samples([4,5,6,7], nkvar)
#export_samples([8,9,10,11,12,13,14,15,16,17,23,24], kt9var)
#export_samples([18,19,20], hhvar)



