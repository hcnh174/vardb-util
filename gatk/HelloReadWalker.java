//package org.vardb.util.tools.nextgen;

import java.io.PrintStream;

import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

/**
 * Walks over the input data set, calculating the total number of covered loci for diagnostic purposes.
 *
 * <p>
 * Simplest example of a locus walker.
 *
 *
 * <h2>Input</h2>
 * <p>
 * One or more BAM files.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * Number of loci traversed.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -cp ./target/vardb-util-0.1.0.BUILD-SNAPSHOT.jar;../lib/GenomeAnalysisTK.jar org.broadinstitute.sting.gatk.CommandLineGATK -T HelloRead -R U:\ref\HCV-KT9.fasta -I U:\OLD2bam\nextgen3-2G__HCV-KT9.bam
 * </pre>
 *
 */
@WalkerName("HelloRead")
public class HelloReadWalker extends ReadWalker<Integer,Integer>
{
	@Output
	public PrintStream out;
	
	@Override
	public Integer map(ReferenceContext ref, GATKSAMRecord read, ReadMetaDataTracker metaDataTracker)
	{
		out.println("Hello, "+read.getReadName());
	    return 1;
	}
	
	@Override
	public Integer reduceInit()
	{
		return null;
	}
	
	@Override
	public Integer reduce(Integer value, Integer sum)
	{
	    return null;
	}
}
