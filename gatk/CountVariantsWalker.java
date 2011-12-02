import java.io.PrintStream;

import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

/**
 * Prints the base counts at each position
 * mvn install && java -Xmx2g -cp ./target/vardb-util-0.1.0.BUILD-SNAPSHOT.jar;../lib/GenomeAnalysisTK.jar org.broadinstitute.sting.gatk.CommandLineGATK -T CountVariants -R U:\ref\HCV-KT9.fasta -I U:\OLD2bam\nextgen3-2G__HCV-KT9.bam -dt NONE -O U:\counts.txt
 * java -Xmx2g -cp ../target/gatk-walkers.jar;../../lib/GenomeAnalysisTK.jar org.broadinstitute.sting.gatk.CommandLineGATK -T CountVariants -R U:\ref\HCV-KT9.fasta -I U:\OLD2bam\nextgen3-2G__HCV-KT9.bam -dt NONE -o U:\counts.txt
 */
//@By(DataSource.REFERENCE)
public class CountVariantsWalker extends LocusWalker<Integer, Integer> implements TreeReducible<Integer>
{
	@Output
	PrintStream out;
	
	public void initialize()
	{
		char tab='\t';
		out.println("position"+tab+"refnt"+tab+"a"+tab+"c"+tab+"g"+tab+"t"+tab+"depth");
	}
	
    @SuppressWarnings("deprecation")
	public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context)
    {
    	ReadBackedPileup basePileup = context.getBasePileup();
    	int[] counts=basePileup.getBaseCounts();
    	out.printf("%d\t",context.getPosition());
    	out.printf("%s\t",ref.getBaseAsChar());
    	out.printf("%d\t%d\t%d\t%d\t",counts[0],counts[1],counts[2],counts[3]);//A, C, G, T
    	out.printf("%d%n",basePileup.depthOfCoverage());  
    	return 1;
    }
    
//    
//    @SuppressWarnings("deprecation")
//	public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context)
//    {
//    	ReadBackedPileup basePileup = context.getBasePileup();
//    	//group	subject	sample	column	ref	region	ntnum	aanum	nt	rank	count	freq
//    	int[] counts=basePileup.getBaseCounts();
//    	out.printf("position=%d",context.getPosition());
//    	out.printf(" ref=%s",ref.getBaseAsChar());
//    	out.printf(" counts=[%d, %d, %d, %d]",counts[0],counts[1],counts[2],counts[3]);
//    	out.printf(" depth=%d",basePileup.depthOfCoverage());
//    	out.printf("%n");
//    	return 1;
//    }

    public Integer reduceInit() { return 0; }
    
    public Integer reduce(Integer value, Integer sum)
    {
        return treeReduce(sum,value);
    }
    
    public Integer treeReduce(Integer lhs, Integer rhs)
    {
        return lhs + rhs;
    }
    
    /*
    // Given result of map function
    @Override
    public void onTraversalDone(Integer result)
    {
        out.println("[REDUCE RESULT] Traversal result is: " + result);
    }
    */
}
