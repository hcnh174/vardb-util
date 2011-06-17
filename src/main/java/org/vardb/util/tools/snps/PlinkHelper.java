package org.vardb.util.tools.snps;

import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Scanner;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.vardb.util.CDataFrame;
import org.vardb.util.CException;
import org.vardb.util.CFileHelper;

public final class PlinkHelper
{
	public static final String SPACE=" ";
	public static final String COMMA=",";
	public static final String TAB="\t";
	public static final String NEWLINE="\n";
	public static final String UNDERSCORE="_";
	
	public static final int BATCH_SIZE=500;

	public static final String FID="FID";
	public static final String IID="IID";
	public static final String MALE="M";
	public static final String FEMALE="F";
	public static final String PATERNAL="0";
	public static final String MATERNAL="0";
	public static final String SR="SR";
	public static final String DISTANCE="0";	
	
	public static final String MAJOR="major";
	public static final String MINOR="minor";
	
	public static final String NO_GENOTYPE="0";
	
	public static final String MISSING="999999999";
	
	public static String stripMissingValueToken(String str)
	{
		return stripMissingValueToken(str,MISSING);
	}
	
	public static String stripMissingValueToken(String str, String missing)
	{
		return str.replace(missing,"");
	}

	public static String convertGenotype(String genotype, String snp_id, CDataFrame snps)
	{
		String major=snps.getStringValue(MAJOR,snp_id);
		String minor=snps.getStringValue(MINOR,snp_id);
		char allele1=genotype.charAt(0);
		char allele2=genotype.charAt(1);
		String str=convertGenotype(allele1,major,minor)+SPACE+convertGenotype(allele2,major,minor);
		//System.out.println("converted genotype "+genotype+" to "+str);
		return str;
	}
	
	private static String convertGenotype(char allele, String major, String minor)
	{
		if (allele=='A')
			return major;
		else if (allele=='B')
			return minor;
		else if (allele=='-')
			return NO_GENOTYPE;
		else { System.out.println("allele is not A or B or -: "+allele); return "0";}
	}
	
	public static String convertGenderPlink(String sex)
	{
		return Sex.find(sex).getPlinkCode().toString();
	}
	
	//0 = female, 1 = male
	public static String convertGenderGenabel(String sex)
	{
		return Sex.find(sex).getGenabelCode().toString();
	}
	
	/*
	public static int convertSexPlink(String sex)
	{
		if (sex==null)
			return 0;
		//System.out.println("sex="+sex);
		if (sex.equals(PlinkHelper.MALE))
			return 1;
		else if (sex.equals(PlinkHelper.FEMALE))
			return 2;
		else return 0; // unknown
	}
	
	//0 = female, 1 = male
	public static String convertSexGenabel(String sex)
	{
		if (sex.equals(PlinkHelper.MALE))
			return "1";
		else if (sex.equals(PlinkHelper.FEMALE))
			return "0";
		throw new CException("can't determine sex for "+sex);
	}
	*/
	
	//Affection status, by default, should be coded:
	//-9 missing 
	//  0 missing
	//1 unaffected
	//2 affected
	public static int convertResponse(String value)
	{
		if (value==null)
			return Integer.parseInt(MISSING);
		return value.equals(PlinkHelper.SR) ? 2 : 1;
	}
	
	public static int convertResponseGenabel(String value)
	{
		return value.equals(PlinkHelper.SR) ? 1 : 0;
	}
	
	public static void writeDataFile(String filename, CDataFrame dataframe)
	{
		String output=dataframe.toString();
		output=PlinkHelper.stripMissingValueToken(output);
		CFileHelper.writeFile(filename,output);
	}
	
	public static void createLDfile(String plinkdir, String bfile, String snp, int window)
	{
		//plink --bfile hcv2 --snp rs8099917 --window 250 --recodeHV --out rs8099917ld250kb"
		System.out.println("creating HaploView LD file for "+snp);
		CCommandLine commands=new CCommandLine("plink");
		commands.setWorkingDir(plinkdir);
		commands.addArg("--bfile",bfile);
		commands.addArg("--snp",snp);
		commands.addArg("--window",window);
		commands.addArg("--recodeHV");
		commands.addArg("--out",snp+"ld"+window+"kb");
		IRuntimeService runtimeService=new RuntimeServiceImpl();
		runtimeService.exec(commands);
	}
	
	public static char flipStrand(char nt)
	{
		if (nt=='A')
			return 'T';
		if (nt=='T')
			return 'A';
		if (nt=='G')
			return 'C';
		if (nt=='C')
			return 'G';
		return nt;
	}
	
	public static Map<String,SnpFrequencyData> getSnpFrequencies(String filename, Collection<String> snpIDs)
	{
		Map<String,SnpFrequencyData> freqs=new LinkedHashMap<String,SnpFrequencyData>();
		for (String snpID : snpIDs)
		{
			SnpFrequencyData freq=new SnpFrequencyData();
			freqs.put(snpID, freq);
		}
		
		Scanner scanner=CFileHelper.createScanner(filename,"\n");
		scanner.nextLine(); //skip the header line
		while(scanner.hasNext())
		{
			String line=scanner.nextLine();
			String[] tabs=line.split("\\s+");
			String snpID=tabs[2];
			if (!freqs.containsKey(snpID))
				continue;
			SnpFrequencyData snp=(SnpFrequencyData)freqs.get(snpID);
			snp.setData(tabs);
		}
		CFileHelper.closeScanner(scanner);
		return freqs;
	}
	
	public static CommandLine parseCommandLine(String[] args, Options options)
	{
		CommandLineParser parser = new BasicParser();//GnuParser();
	    try
	    {
	        CommandLine line = parser.parse(options, args);
	        return line;
	    }
	    catch(ParseException e)
	    {
	       throw new CException(e);
	    }
	}
	
	public static class SnpFrequencyData
	{
		public String chr;
		public String snp;
		public String a1;
		public String a2;
		public String maf;
		public String nchrobs;
		
		public void setData(String[] tabs)
		{
			chr=tabs[1];
			snp=tabs[2];
			a1=tabs[3];
			a2=tabs[4];
			maf=tabs[5];
			nchrobs=tabs[6];
			//System.out.println("setting freq: SNP="+snp+" maf="+maf+" "+nchrobs);
		}
	}
	
	public enum Sex
	{
		MALE("M",1,1),
		FEMALE("F",2,0);
		
		//plink  sex (1=male; 2=female; other=unknown)
		//genabel (0 = female, 1 = male)
		
		Sex(String letter, int plink, int genabel)
		{
			this.letter=letter;
			this.plink=plink;
			this.genabel=genabel;
		}
		
		private String letter;
		private int plink;
		private int genabel;
		
		public String getLetter(){return this.letter;}
		public Integer getPlinkCode(){return this.plink;}
		public Integer getGenabelCode(){return this.genabel;}
		
		public static Sex find(String value)
		{
			value=value.trim();
			for (Sex sex : values())
			{
				if (sex.getLetter().equals(value))// || sex.getCode().toString().equals(value))
					return sex;
			}
			throw new CException("Can't determine sex for value: "+value);
		}
	}
	
	/*
	-9 missing 
    0 missing
    1 unaffected
    2 affected
	*/
	public enum Response
	{
		AFFECTED("2","1"),
		UNAFFECTED("1","0"),
		MISSING("0","NA");
		
		Response(String plink, String genabel)
		{
			this.plink=plink;
			this.genabel=genabel;
		}
		private String plink;
		private String genabel;
		
		public String getPlinkCode(){return this.plink;}
		public String getGenabelCode(){return this.genabel;}
		
		public static Response findByPlinkCode(String value)
		{
			value=value.trim();
			if (value.equals(PlinkHelper.MISSING))
				return MISSING;
			for (Response response : values())
			{
				if (response.getPlinkCode().equals(value))
					return response;
			}
			throw new CException("Can't determine response for value: "+value);
		}
	}
}