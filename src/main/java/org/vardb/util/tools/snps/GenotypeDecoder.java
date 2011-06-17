package org.vardb.util.tools.snps;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.vardb.util.CDataFrame;
import org.vardb.util.CException;
import org.vardb.util.CFileHelper;
import org.vardb.util.CStringHelper;

public final class GenotypeDecoder
{
	public static void main(String[] args)
	{
		/*
		Database db=new Database();
		
		ExpressionParser parser=new SpelExpressionParser();
		Expression exp=parser.parseExpression("name");
		EvaluationContext context=new StandardEvaluationContext(db);
		//((StandardEvaluationContext) context).setRootObject(db);
		System.out.println("name="+(String)exp.getValue(context));
		
		System.out.println("db="+(String)exp.getValue(db));
		*/
		
		Params params=new Params(args);
		System.out.println("params="+params.toString());
		
		GenotypeDecoder decoder=new GenotypeDecoder(params);
		CDataFrame dataframe=decoder.decode();
		CFileHelper.writeFile(params.outfile,dataframe.toString());
	}
	
	/*
	public static class Database
	{
		public String getName(){return "Nelson";}
	}
	*/
	/////////////////////////////////////////////////////////////////////////
	
	protected Params params;
	
	public GenotypeDecoder(Params params)
	{
		this.params=params;
	}
	
	public CDataFrame decode()
	{
		CDataFrame alleles=loadAlleleFile(params.allelefile);
		CDataFrame data=CDataFrame.parseTabFile(params.datafile);
		List<String> snpnames=getSnpNames(data);
		for (String snpname : snpnames)
		{
			CDataFrame.Column column=data.getColumn(snpname);
			String a=alleles.getStringValue("A",snpname);
			String b=alleles.getStringValue("B",snpname);
			for (Object rowname : column.getRowNames())
			{
				String value=column.getValue(rowname).toString();
				String newvalue=value;
				if (value.equals("allele1") || value.equals("AA"))
					newvalue=a+a;
				else if (value.equals("Both") || value.equals("AB"))
					newvalue=a+b;
				else if (value.equals("allele2") || value.equals("BB"))
					newvalue=b+b;
				else if (value.equals("Undetermined"))
					newvalue="";
				column.setValue(rowname, newvalue);
			}
		}
		return data;
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
		public String outfile="data-decoded.txt";
		
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
			this.outfile=this.datafile.replace(".txt", "-decoded.txt");
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
