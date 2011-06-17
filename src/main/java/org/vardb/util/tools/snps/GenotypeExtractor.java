package org.vardb.util.tools.snps;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;

import org.vardb.util.CDataFrame;
import org.vardb.util.CException;
import org.vardb.util.CFileHelper;
import org.vardb.util.CStringHelper;

// given a set of SNP IDs, creates a table with the patient ID in the first column and then one column for each SNP
public final class GenotypeExtractor extends CAbstractCommandLineWrapper
{
	
	public static void main(String[] args)
	{
		String plinkdir=CFileHelper.normalizeDirectory(args[0]);
		String bfile=args[1];
		String phenofile=args[2];
		String outfile=args[3];
		List<String> snps=CStringHelper.splitAsList(args[4],",");

		boolean indicators=true;		
		CDataFrame dataframe=CDataFrame.parseTabFile(plinkdir+phenofile);

		GenotypeExtractor extractor=new GenotypeExtractor(new Params(plinkdir,bfile));
		CDataFrame genotypes=extractor.extract(snps,indicators);			
		dataframe.appendColumns(genotypes);
	
		PlinkHelper.writeDataFile(plinkdir+outfile,dataframe);
	}
	
	private Params params;
		
	public GenotypeExtractor(Params params)
	{
		super(params.dir,params.dir);
		this.params=params;
	}
	
	public CDataFrame extract(Collection<String> snpIDs, boolean indicators)
	{
		if (snpIDs.isEmpty())
			throw new CException("no SNPs to look up");
		System.out.println("extracting SNPs: "+CStringHelper.join(snpIDs,", "));
		Params params=new Params(this.params.dir,this.params.bfile);
		String prefix="snplist";
		params.out=getRoot(prefix);
		params.infile=createTempfile(prefix,".list",CStringHelper.join(snpIDs,"\n"));
		System.out.println("extracting SNP genotypes from directory "+params.dir+" in list "+params.infile);		
		addTempfile(params.dir+params.out+".log");
		addTempfile(params.dir+params.out+".map");
		addTempfile(params.dir+params.out+".ped");
		
		CCommandLine commands=getCommand(params);
		runtimeService.exec(commands);
		
		List<Snp> snps=getSnps(params);
		List<Genotype> genotypes=getGenotypes(params,snps);
		CDataFrame dataframe=createReport(snps,genotypes,indicators);
		deleteTempfiles();
		return dataframe;
	}
	
	public CDataFrame extractFromFile(String filename, boolean indicators)
	{
		if (!CFileHelper.exists(params.dir+params.infile))
			throw new CException("can't find extract file: "+params.dir+params.infile);
		List<String> snpIDs=CStringHelper.splitLines(CFileHelper.readFile(filename));
		return extract(snpIDs,indicators);
	}
	
	private static CCommandLine getCommand(Params params)
	{
		//String cmd=params.dir+"plink --bfile hcv --recode --extract "+params.listfilename+" --out genotypes";
		CCommandLine commands=new CCommandLine("plink");
		commands.setWorkingDir(params.dir);
		commands.addArg("--bfile",params.bfile);
		commands.addArg("--recode");
		commands.addArg("--extract",params.infile);
		commands.addArg("--out",params.out);
		return commands;
	}
	
	private static CDataFrame createReport(List<Snp> snps, List<Genotype> genotypes, boolean indicators)
	{
		CDataFrame dataframe=new CDataFrame();
		dataframe.addColumn(PlinkHelper.IID);
		for (Snp snp : snps)
		{
			dataframe.addColumn(snp.snp);
			if (indicators)
			{
				dataframe.addColumn(snp.snp+"ABC");
				dataframe.addColumn(snp.snp+"AA");
				dataframe.addColumn(snp.snp+"AB");
				dataframe.addColumn(snp.snp+"BB");
				dataframe.addColumn(snp.snp+"Adom");
				dataframe.addColumn(snp.snp+"Bdom");
			}
		}
		for (Genotype genotype : genotypes)
		{
			dataframe.setValue(PlinkHelper.IID,genotype.iid,genotype.iid);
			for (Snp snp : snps)
			{
				setGenotype(dataframe,snp.snp,genotype.iid,genotype.genotypes.get(snp),indicators);
			}
		}
		return dataframe;
	}
	
	private static void setGenotype(CDataFrame dataframe, String snp, String iid, String genotype, boolean indicator)
	{
		if (genotype.equals("00"))
			genotype="";
		dataframe.setValue(snp,iid,genotype);
		if (!indicator)
			return;
		int abc=0;
		int aa=0;
		int ab=0;
		int bb=0;
		int adom=0;
		int bdom=0;
		if (genotype.equals("AA"))
		{
			aa=1;
			adom=1;
			abc=1;
		}
		else if (genotype.equals("AB") || genotype.equals("BA"))
		{
			ab=1;
			adom=1;
			bdom=1;
			abc=2;
		}
		else if (genotype.equals("BB"))
		{
			bb=1;
			bdom=1;
			abc=3;
		}
		
		dataframe.setValue(snp+"ABC",iid,abc);
		dataframe.setValue(snp+"AA",iid,aa);
		dataframe.setValue(snp+"AB",iid,ab);
		dataframe.setValue(snp+"BB",iid,bb);
		dataframe.setValue(snp+"Adom",iid,adom);
		dataframe.setValue(snp+"Bdom",iid,bdom);
	}
	
	//7	rs1523638	0	152267562
	private static List<Snp> getSnps(Params params)
	{
		List<Snp> snps=new ArrayList<Snp>();
		Scanner scanner=CFileHelper.createScanner(params.dir+params.out+".map");
		while (scanner.hasNext())
		{
			String line=scanner.next();
			Snp snp=new Snp(line);
			snps.add(snp);
		}
		CFileHelper.closeScanner(scanner);
		return snps;
	}
	
	//DD01_A01 DD01_A01 0 0 1 1 B B A A
	private static List<Genotype> getGenotypes(Params params, List<Snp> snps)
	{
		List<Genotype> genotypes=new ArrayList<Genotype>();
		Scanner scanner=CFileHelper.createScanner(params.dir+params.out+".ped");
		while (scanner.hasNext())
		{
			String line=scanner.next();
			Genotype genotype=new Genotype(line,snps);
			genotypes.add(genotype);
		}
		CFileHelper.closeScanner(scanner);
		return genotypes;
	}
	
	public static class Genotype
	{
		public String iid;
		public Map<Snp,String> genotypes=new LinkedHashMap<Snp,String>();
		
		public Genotype(String line, List<Snp> snps)
		{
			line=line.trim();
			String[] tabs=line.split("\\s+");
			//System.out.println("line=["+line+"]");
			iid=tabs[1];			
			int snpindex=0;
			for (int index=6;index<tabs.length;index+=2)
			{
				Snp snp=snps.get(snpindex);
				genotypes.put(snp,tabs[index]+tabs[index+1]);
				snpindex++;
			}
		}
		
		public String getGenotype(Snp snp)
		{
			return genotypes.get(snp);
		}
	}
	
	public static class Snp
	{
		public String chr;
		public String snp;
		public String bp;
		
		public Snp(String line)
		{
			line=line.trim();
			//System.out.println("line=["+line+"]");
			String[] tabs=line.split("\\s+");
			chr=tabs[0];
			snp=tabs[1];
			bp=tabs[3];
		}
	}
	
	public static class Params
	{
		public String dir;
		public String bfile;
		public String infile;
		public String out="snplist";
		public String outfile="genotypes.geno";
		public String phenotypefile="pheno.phe";
		public String genotypefile="genotypes.geno";
		public String phengenoutfile=phenotypefile+"+geno.txt";
		public String missing="999999999";
		
		public Params(String dir, String bfile)
		{
			this.dir=CFileHelper.normalizeDirectory(dir);
			this.bfile=bfile;
		}
	}
}