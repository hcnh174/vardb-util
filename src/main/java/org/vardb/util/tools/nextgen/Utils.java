package org.vardb.util.tools.nextgen;

import org.vardb.util.CTable;

//read in directories and sample information
// convert all files to fasta or at least fastq
// trim reads using quality control
// run velveth along with reference sequences
// run velvetg
// align contigs with reference sequence
// call snps
public class Utils
{
	public static void main(String args[])
	{
		//loadSampleSheet(args[0]);
	}

	public SampleSheet loadSampleSheet(String filename)
	{
		SampleSheet samplesheet=new SampleSheet();
		CTable.parseFile(filename,samplesheet,",");
		return samplesheet;
	}
	
	public static class SampleSheet extends CTable
	{
		
	}
}
