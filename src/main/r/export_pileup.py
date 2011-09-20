import sys
import pysam

def export_sample(sample, ref, bamdir, pileupdir):
	print 'sample='+sample+', ref='+ref+', bamdir='+bamdir+', pileupdir='+pileupdir
	#infile = "bam/"+sample+"."+ref+".bam"
	infile =bamdir+"/"+sample+".bam"
	outfile = pileupdir+"/"+sample+".txt" 
	print 'infile='+infile+', outfile='+outfile
	out = open(outfile, "w")
	out.write('position\tread\tnt\n')
	samfile = pysam.Samfile(infile, "rb" )
	for pileupcolumn in samfile.pileup(ref):
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

# python export_pileup.py sample18 HHaa36
sample = sys.argv[1] 
ref = sys.argv[2]
bamdir = sys.argv[3]
pileupdir = sys.argv[4]

export_sample(sample,ref,bamdir,pileupdir)
