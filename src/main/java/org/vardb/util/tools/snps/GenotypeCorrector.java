package org.vardb.util.tools.snps;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.vardb.util.CDataFrame;
import org.vardb.util.CException;
import org.vardb.util.CFileHelper;
import org.vardb.util.CStringHelper;

import com.google.common.collect.Maps;

public final class GenotypeCorrector
{
	public static void main(String[] args)
	{
		Params params=new Params(args);
		System.out.println("params="+params.toString());
		
		GenotypeCorrector corrector=new GenotypeCorrector(params);
		CDataFrame dataframe=corrector.correct();
		CFileHelper.writeFile(params.outfile,dataframe.toString());
	}
	
	/////////////////////////////////////////////////////////////////////////
	
	protected Params params;
	private Map<Character,Character> bases=Maps.newHashMap();
	
	public GenotypeCorrector(Params params)
	{
		this.params=params;
		bases.put('A','T');
		bases.put('T','A');
		bases.put('C','G');
		bases.put('G','C');
	}

	public CDataFrame correct()
	{
		CDataFrame alleles=loadAlleleFile(params.allelefile);
		CDataFrame data=CDataFrame.parseTabFile(params.datafile);
		List<String> snpnames=getSnpNames(data);
		for (String snpname : snpnames)
		{
			CDataFrame.Column column=data.getColumn(snpname);
			char a=alleles.getStringValue("A",snpname).charAt(0);
			char b=alleles.getStringValue("B",snpname).charAt(0);
			for (Object rowname : column.getRowNames())
			{
				String value=column.getValue(rowname).toString();
				if (value.equals("NA") || value.equals(""))
					continue;
				char allele1=value.charAt(0);
				char allele2=value.charAt(1);
				if (alleleOkay(allele1,a,b) & alleleOkay(allele2,a,b))
					column.setValue(rowname,makeGenotype(allele1,allele2));
				else column.setValue(rowname,flipGenotype(allele1,allele2));
			}
		}
		return data;
	}
	
	private String flipGenotype(char allele1, char allele2)
	{
		allele1=flipStrand(allele1);
		allele2=flipStrand(allele2);
		return makeGenotype(allele1,allele2);
	}
	
	private String makeGenotype(char allele1, char allele2)
	{
		if (allele1>allele2)
			return ""+allele2+""+allele1;
		else return ""+allele1+""+allele2;
	}
		
	private char flipStrand(char allele)
	{
		return bases.get(allele);
	}
	
	private boolean alleleOkay(char allele, char a, char b)
	{
		return (allele==a || allele==b); 
	}
	
	private List<String> getSnpNames(CDataFrame data)
	{
		List<String> snpnames=new ArrayList<String>();
		for (String colname : data.getColNames())
		{
			if (colname.indexOf("rs")==0)
				snpnames.add(colname);
		}
		return snpnames;
	}
	
	private CDataFrame loadAlleleFile(String filename)
	{
		CDataFrame alleles=CDataFrame.parseTabFile(filename);
		alleles.addColumn("A");
		alleles.addColumn("B");
		CDataFrame.Column column=alleles.getColumn("alleles");
		for (Object rowname : column.getRowNames())
		{
			String value=column.getValue(rowname).toString();
			String a=""+value.charAt(0);
			String b=""+value.charAt(value.length()-1);
			alleles.setValue("A",rowname,a);
			alleles.setValue("B",rowname,b);
		}
		return alleles;
	}
	
	public static class Params
	{
		public String dir="c:\\temp\\";
		public String datafile="data.txt";
		public String allelefile="alleles.txt";
		public String outfile="data-corrected.txt";
		
		public Params(){}
		
		public Params(String[] args)
		{
			setParams(args);
		}
		
		public void setParams(String[] argv)
		{
			Options options = new Options();
			options.addOption("dir", true,"input directory");
			options.addOption("datafile", true,"data input file");
			options.addOption("allelefile", true,"allele lookup table");
			//options.addOption("outfile", true,"output file");
			
			CommandLine line=PlinkHelper.parseCommandLine(argv,options);
			if (line.hasOption("dir")) this.dir=CFileHelper.normalizeDirectory(line.getOptionValue("dir"));
			if (line.hasOption("datafile")) this.datafile=fixPath(line.getOptionValue("datafile"));
			if (line.hasOption("allelefile")) this.allelefile=fixPath(line.getOptionValue("allelefile"));
			//if (line.hasOption("outfile")) this.outfile=fixPath(line.getOptionValue("outfile"));
			//if (line.hasOption("outfile")) this.outfile=line.getOptionValue("outfile");
			this.outfile=this.datafile.replace(".txt", "-corrected.txt");
		}
		
		private String fixPath(String filename)
		{
			filename=CFileHelper.normalize(filename);
			if (filename.charAt(0)!='/' && filename.charAt(1)!=':') // either unix-style or c: style
				filename=this.dir+filename;
			if (!CFileHelper.exists(filename))
				throw new CException("can't find filename: "+filename);
			return filename;
		}
		
		@Override
		public String toString()
		{
			return CStringHelper.toString(this);
		}
		

	}
}
