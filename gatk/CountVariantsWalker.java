import java.io.PrintStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.IntervalBinding;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.arguments.GATKArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

/**
 * Prints the base counts at each position
 * mvn install && java -Xmx2g -cp ./target/vardb-util-0.1.0.BUILD-SNAPSHOT.jar;../lib/GenomeAnalysisTK.jar org.broadinstitute.sting.gatk.CommandLineGATK -T CountVariants -R U:\ref\HCV-KT9.fasta -I U:\OLD2bam\nextgen3-2G__HCV-KT9.bam -dt NONE -O U:\counts.txt
 * java -Xmx10g -cp ../target/gatk-walkers.jar;../../lib/GenomeAnalysisTK.jar org.broadinstitute.sting.gatk.CommandLineGATK -T CountVariants -R U:\ref\HCV-KT9.fasta -I U:\OLD2bam\nextgen3-2G__HCV-KT9.bam -dt NONE -o U:\counts.txt -L HCV-KT9:6139-6174
 */
//@By(DataSource.REFERENCE)
public class CountVariantsWalker extends LocusWalker<Integer, Integer> implements TreeReducible<Integer>
{
	@Output
	protected PrintStream out;
	
	protected Integer baseQuality=30;
	protected Integer mapQuality=30;
	
	protected Integer pos=0;
	protected Map<String,Character> nt1=new HashMap<String,Character>();
	protected Map<String,Character> nt2=new HashMap<String,Character>();
	protected Map<String,Character> nt3=new HashMap<String,Character>();
	
	public void initialize()
	{
		char tab='\t';
		out.println("position"+tab+"refnt"+tab+"a"+tab+"c"+tab+"g"+tab+"t"+tab+"depth");
		GATKArgumentCollection args=super.getToolkit().getArguments();
//    	for (IntervalBinding<Feature> interval : args.intervals)
//    	{
//    		System.out.println("interval: "+interval.toString());
//    	}
	}
	
	/*
	private boolean isNewInterval(int position)
	{
		GATKArgumentCollection args=super.getToolkit().getArguments();
		for (IntervalBinding<Feature> interval : args.intervals)
		{
			
		}
	}
	*/
	
	private Map<String,Character> getNtMap()
	{
    	if (pos%3==0)
    		return nt1;
    	else if (pos%3==1)
    		return nt2;
    	else if (pos%3==2)
    		return nt3;
    	else throw new RuntimeException("unexpected position: "+pos);
	}
	
    //@SuppressWarnings("deprecation")
	public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context)
    {
    	System.out.println("position="+pos);
    	Map<String,Character> nt=getNtMap();
    	
    	//System.out.println("*********************************************");
    	//System.out.println("position="+context.getPosition());
    	//ReadBackedPileup basePileup = context.getBasePileup();
    	ReadBackedPileup pileup = context.getBasePileup();
    	ReadBackedPileup basePileup = pileup.getBaseAndMappingFilteredPileup(baseQuality,mapQuality);
    	List<GATKSAMRecord> reads=basePileup.getReads();
    	byte[] bases=basePileup.getBases();
    	if (reads.size() != bases.length)
    		throw new RuntimeException("reads.size="+reads.size()+" is not equal to bases.length="+bases.length);
    	//System.out.println("reads.size="+reads.size()+", bases.length="+bases.length);
    	//for (GATKSAMRecord read : reads)
    	for (int index=0; index<reads.size();index++)
    	{
    		GATKSAMRecord read=reads.get(index);
    		char base=(char)bases[index];
    		nt.put(read.getReadName(),base);
    	}
    	//ReadBackedPileup getBaseAndMappingFilteredPileup( int minBaseQ, int minMapQ )
    	//List<GATKSAMRecord> getReads();
    	//List<Integer> getOffsets();
    	//ReadBackedPileup getFilteredPileup(PileupElementFilter filter)
    	int[] counts=basePileup.getBaseCounts();
    	out.printf("%d\t",context.getPosition());
    	out.printf("%s\t",(char)ref.getBase());
    	out.printf("%d\t%d\t%d\t%d\t",counts[0],counts[1],counts[2],counts[3]);//A, C, G, T
    	out.printf("%d%n",basePileup.depthOfCoverage()); 
    	
    	if (pos%3==2)
    	{
    		Map<String,Integer> codons=new HashMap<String,Integer>();
    		for (String name : nt1.keySet())
    		{
    			if (nt2.containsKey(name) && nt3.containsKey(name))
    			{
    				String codon=""+nt1.get(name)+nt2.get(name)+nt3.get(name);
    				//System.out.println("codon="+codon);
    				if (!codons.containsKey(codon))
    					codons.put(codon,1);
    				else codons.put(codon,codons.get(codon)+1);
    			}
    		}
    		System.out.println("*********************************************");
    		for (String codon : codons.keySet())
    		{
    			System.out.println("codon="+codon+"="+codons.get(codon));
    		}
    		nt1.clear();
    		nt2.clear();
    		nt3.clear();
    	}
    	pos++;
    	return 1;
    }

    public Integer reduceInit() { return 0; }
    
    public Integer reduce(Integer value, Integer sum)
    {
        return treeReduce(sum,value);
    }
    
    public Integer treeReduce(Integer lhs, Integer rhs)
    {
        return lhs + rhs;
    }
}
