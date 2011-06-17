package org.vardb.util.tools.snps;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Scanner;
import java.util.Set;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.vardb.util.CDataFrame;
import org.vardb.util.CException;
import org.vardb.util.CFileHelper;
import org.vardb.util.CStringHelper;
import org.vardb.util.CTable;
import org.vardb.util.tools.snps.PlinkHelper.Response;

import com.google.common.io.Files;
import com.google.common.io.LineProcessor;

public final class PlinkConverter
{
	private static final String COLUMN_FID="FID";
	private static final String COLUMN_IID="IID";
	private static final String COLUMN_GENDER="gender";//"sex";//"Gender";
	private static final String COLUMN_SEX="sex";//"sex";//"Gender";
	private static final String COLUMN_AGE="age";//"Age";
	private static final String COLUMN_CLASS="class";
	private static final String COLUMN_TREATMENT="treatment";
	private static final String COLUMN_NUMTREATMENTS="numtreatments";
	private static final String COLUMN_RESPONSE="response";

	private final static String BOT="BOT";
	private final static String NA="NA";
	private final static boolean shortcut=false;
	private final static int GENOTYPE_START_INDEX=4;
		
	public static void main(String[] args)
	{
		Params params=new Params(args);
		System.out.println("params="+params.toString());
		
		//CPlinkConverter2.convertSampleFile(params.samplefile, params.classfile, "c:/temp/pheno.phe");
		//if (true)
		//	return;
		
		PlinkConverter converter=new PlinkConverter(params);
		converter.getPatientIdentifiers();
		converter.createGenabelPhenofile();

		//if (shortcut) return;
		converter.createTpedFile();
		converter.createTfamFile();
		converter.createBfile();
	}
	
	/////////////////////////////////////////////////////////////////////////
	
	protected Params params;
	protected List<String> patient_ids;
	protected CDataFrame phenotypes;
	protected Set<String> snplist;
	
	public PlinkConverter(Params params)
	{
		this.params=params;
		getPhenotypeDataFrame();
	}
	
	private void loadSnplist()
	{
		if (params.snplistfile==null)
			return;
		snplist=new HashSet<String>();
		try		
		{
			Files.readLines(new File(params.snplistfile), CFileHelper.ENCODING, new LineProcessor<String>()
			{
				public boolean processLine(String line)
				{
					snplist.add(line);
					return true;
				}
				
				public String getResult(){return null;};
			});	
		}
		catch (IOException e)
		{
			throw new CException(e);
		}
	}
	
	//family ID
	//individual ID
	//snp ID
	//allele 1 of this genotype
	//allele 2 of this genotype
	public void createTpedFile()
	{
		String outfile=params.outdir+params.name+".tped";
		System.out.println("creating tped file from "+params.snpfile+" and "+params.genofile+": "+outfile);
		loadSnplist();
		CDataFrame snps=parseSnpFile(params.snpfile);
		Scanner scanner=CFileHelper.createScanner(params.genofile,"\n");
		moveToData(scanner);		
		//this.patient_ids=getPatientIdentifiers(scanner);
		if (shortcut)return;
		int counter=0;
		StringBuilder buffer=new StringBuilder();
		CFileHelper.writeFile(outfile);
		while(scanner.hasNext())
		{
			String line=scanner.nextLine();
			if (line.startsWith("SNP,Chr,Position,Allele"))
				continue;
			if (line.trim().equals(""))
				continue;
			//System.out.println("line: "+line);
			List<String> tabs=CStringHelper.splitAsList(line,PlinkHelper.COMMA);
			//System.out.println("tabs.size="+tabs.size()); 
			if (!(tabs.size()==GENOTYPE_START_INDEX+patient_ids.size()))
			{
				System.out.println("tabs.size is does not match expected size ("+tabs.size()+" vs "+GENOTYPE_START_INDEX+patient_ids.size()+"): ["+line+"]");
				continue;
			}
			
			String snp_id=tabs.get(0);
			if (!keepSnp(snp_id))
			{
				//System.out.println("skipping snp: "+snp_id);
				continue;
			}
			//List<String> genotypes=tabs.subList(1,tabs.size());
			List<String> genotypes=tabs.subList(GENOTYPE_START_INDEX,tabs.size());
			if (genotypes.size()!=patient_ids.size())
				throw new CException("Genotypes.size != patient_ids.size ("+genotypes.size()+" vs "+patient_ids.size());
			
			//System.out.println("chr.tabs="+tabs.get(1)+" vs chr.lookup="+snps.getStringValue("Chr",snp_id));
			
			buffer.append(snps.getStringValue("Chr",snp_id)).append(PlinkHelper.TAB);
			buffer.append(snp_id).append(PlinkHelper.TAB);
			buffer.append(PlinkHelper.DISTANCE).append(PlinkHelper.TAB);
			buffer.append(snps.getStringValue("MapInfo",snp_id));
			
			for (int index=0; index<genotypes.size(); index++)
			{
				String genotype=genotypes.get(index);
				String patient_id=patient_ids.get(index);
				if (!keepSample(patient_id))
				{
					//System.out.println("skipping patient "+patient_id+", index="+index);
					continue;
				}
				buffer.append(PlinkHelper.TAB);
				buffer.append(PlinkHelper.convertGenotype(genotype,snp_id,snps));
			}
			buffer.append(PlinkHelper.NEWLINE);
			
			if (counter%PlinkHelper.BATCH_SIZE==0)
			{
				if (counter%10000==0)
					System.out.println("writing buffer="+counter+" rows");
				CFileHelper.appendFile(outfile,buffer.toString(),false);
				buffer=new StringBuilder();
			}
			counter++;
		}
	}
	
	private boolean keepSnp(String snp_id)
	{
		if (params.rsonly && snp_id.indexOf("rs")==-1)
			return false;
		if (params.snplistfile!=null && !snplist.contains(snp_id))
		{
			//System.out.println("snp not found in snplist, skipping: "+snp_id);
			return false;
		}
		return true;
	}

	private boolean keepSample(Object sample_id)
	{
		if (params.minimize && !phenotypes.hasRow(sample_id))
		//if (!phenotypes.hasRow(sample_id))
			return false;
		return true;
	}
	
	private static CDataFrame parseSnpFile(String snpfile)
	{
		System.out.println("parsing SNP file: "+snpfile);
		CDataFrame.TabFileParser parser=new CDataFrame.TabFileParser()
		{
			@Override
			protected List<String> preProcessHeader(List<String> fields)
			{
				System.out.println("pre-processing header: "+CStringHelper.join(fields,"|"));
				fields.add(PlinkHelper.MAJOR);
				fields.add(PlinkHelper.MINOR);
				return fields;
			}
			
			@Override
			protected List<String> preProcessLine(List<String> values, int rownum)
			{
				String alleles=values.get(1); //hack! - try not to hard-code indexes
				String strand=values.get(4);
				char major=alleles.charAt(1);
				char minor=alleles.charAt(3);
				if (BOT.equals(strand))
				{
					major=PlinkHelper.flipStrand(major);
					minor=PlinkHelper.flipStrand(minor);
				}
				values.add(String.valueOf(major));
				values.add(String.valueOf(minor));
				//System.out.println("pre-processed line "+rownum+": "+rowname+": "+CStringHelper.join(values,"|"));
				return values;
			}
		};
		parser.parseFile(snpfile);
		CDataFrame snps=parser.getDataFrame();
		System.out.println("finished parsing SNP file");
		//CFileHelper.writeFile("c:/temp/alleles.txt",snps.toString());		
		return snps;
	}
	
	private static void moveToData(Scanner scanner)
	{
		int counter=0;
		while (true)
		{
			String line=scanner.nextLine();
			counter++;
			if (counter>15)
				throw new CException("already past 15 lines but did not find [DATA]. last line="+line);
			if ("[Data]".equals(line))
				return;
		}
	}

	/*
	private static List<String> getPatientIdentifiers(Scanner scanner)
	{
		System.out.println("getting patient identifiers");
		String line=scanner.nextLine();
		List<String> ids=CStringHelper.splitAsList(line,PlinkHelper.COMMA);
		ids=ids.subList(1,ids.size());
		System.out.println("done getting patient identifiers");
		return ids;
	}
	*/
	
	public void getPatientIdentifiers()
	{
		System.out.println("getting patient identifiers");
		Scanner scanner=CFileHelper.createScanner(params.genofile,"\n");
		moveToData(scanner);		
		String line=scanner.nextLine();
		List<String> ids=CStringHelper.splitAsList(line,PlinkHelper.COMMA);
		ids=ids.subList(GENOTYPE_START_INDEX,ids.size());
		System.out.println("done getting patient identifiers");
		this.patient_ids=ids;
	}
	
	////////////////////////////////////////////////////////////////////
	
	//Family ID
    //Individual ID
    //Paternal ID
    //Maternal ID
    //Sex (1=male; 2=female; other=unknown)
    //Phenotype
	public void createTfamFile()
	{
		System.out.println("creating tfam file from phenofile: "+params.phenofile);
		CTable table=new CTable();
		table.setShowHeader(false);
		table.getHeader().add("Family ID");
		table.getHeader().add("Individual ID");
		table.getHeader().add("Paternal ID");
		table.getHeader().add("Maternal ID");
		table.getHeader().add("Sex");
		table.getHeader().add("Phenotype");

		for (Object patient_id : patient_ids)
		{
			//System.out.println("patient "+patient_id);
			//if (!phenotypes.hasRow(patient_id))
			//{
				//System.out.println("no data for patient: "+patient_id);
			//	continue;
			//}
			if (!keepSample(patient_id))
			{
				//System.out.println("skipping patient "+patient_id);
				continue;
			}
			CTable.Row row=table.addRow();
			row.add(patient_id);
			row.add(patient_id);
			row.add(PlinkHelper.PATERNAL);
			row.add(PlinkHelper.MATERNAL);
			row.add(PlinkHelper.convertGenderPlink(phenotypes.getStringValue(COLUMN_GENDER,patient_id)));
			row.add(PlinkHelper.convertResponse(phenotypes.getStringValue(COLUMN_RESPONSE,patient_id)));
		}
		String outfile=params.outdir+params.name+".tfam";
		CFileHelper.writeFile(outfile,table.toString());
	}
	
	public void createGenabelPhenofile()
	{
		CDataFrame dataframe=CDataFrame.parseTabFile(params.phenofile);
		CDataFrame.Column gendercolumn=dataframe.getColumn(COLUMN_GENDER);
		CDataFrame.Column sexcolumn=dataframe.getColumn(COLUMN_SEX);
		for (Object rowname : gendercolumn.getRowNames())
		{
			String gender=gendercolumn.getValue(rowname).toString();
			sexcolumn.setValue(rowname, PlinkHelper.convertGenderGenabel(gender));
		}
		CDataFrame.Column column=dataframe.getColumn("response");
		CDataFrame.Column svrcolumn=dataframe.addColumn("svr");
		CDataFrame.Column nrcolumn=dataframe.addColumn("nr");
		for (Object rowname : column.getRowNames())
		{
			String response=column.getValue(rowname).toString();
			svrcolumn.setValue(rowname, makeSvrField(response));
			nrcolumn.setValue(rowname, makeNrField(response));
		}		
		String str=dataframe.toString();
		str=str.replace(PlinkHelper.MISSING, NA);
		str=str.replace(PlinkHelper.FID, "id");
		String filename=params.outdir+params.name+".pheno";
		System.out.println("writing phenofile: "+filename);
		CFileHelper.writeFile(filename,str);
	}
	
	private String makeSvrField(String value)
	{
		Response response=Response.findByPlinkCode(value);
		if (response==Response.AFFECTED)
			return "0";
		else if (response==Response.UNAFFECTED)
			return "1";
		else return PlinkHelper.MISSING;
	}
	
	private String makeNrField(String value)
	{
		Response response=Response.findByPlinkCode(value);
		return response.getGenabelCode();
	}
	
	/*
	private String convertSvrGenabel(String response)
	{
		if (response.equals("SVR"))
			return "1";
		else if (response.equals("TR") || response.equals("NR"))
			return "0";
		else return "NA";
	}
	
	private String convertNrGenabel(String response)
	{
		if (response.equals("NR"))
			return "1";
		else if (response.equals("TR") || response.equals("SVR"))
			return "0";
		else return "NA";
	}
	 */
	
	private CDataFrame getPhenotypeDataFrame()
	{
		if (phenotypes==null)
			phenotypes=CDataFrame.parseTabFile(params.phenofile);
		return phenotypes;
	}
	
	/*
	ID      96plate loc     class   QC check        Gender  Age
	H431    DD01    A01     <4>             M       64
	H505    DD01    A02     <7>             M       82
	*/	
	public static void convertSampleFile(String samplefile, String classfile, String outfile)
	{	
		CDataFrame classes=CDataFrame.parseTabFile(classfile);
		CDataFrame dataframe=new CDataFrame();
		dataframe.addColumn(COLUMN_FID);
		dataframe.addColumn(COLUMN_IID);
		dataframe.addColumn(COLUMN_CLASS);
		dataframe.addColumn(COLUMN_SEX);
		dataframe.addColumn(COLUMN_AGE);
		dataframe.addColumn(COLUMN_TREATMENT);
		dataframe.addColumn(COLUMN_NUMTREATMENTS);
		dataframe.addColumn(COLUMN_RESPONSE);
		Scanner scanner=CFileHelper.createScanner(samplefile,"\n");
		// skip first line
		scanner.nextLine();
		while(scanner.hasNext())
		{
			String line=scanner.nextLine();
			List<String> tabs=CStringHelper.splitAsList(line,PlinkHelper.TAB);
			String patient_id=tabs.get(1)+PlinkHelper.UNDERSCORE+tabs.get(2);
			
			String cls=tabs.get(3);
			cls=cls.substring(1,cls.length()-1);
			if (!classes.hasValue("treatment",cls))
				throw new CException("cannot find treatment value for class "+cls);
			
			String sex=tabs.get(5);
			//int sex=PlinkHelper.convertSex(tabs.get(5));
			String age=tabs.get(6);
			String treatment=classes.getStringValue(COLUMN_TREATMENT,cls);
			int numtreatments=classes.getIntValue(COLUMN_NUMTREATMENTS,cls);
			String response=classes.getStringValue(COLUMN_RESPONSE,cls);
			//int response=PlinkHelper.convertResponse(classes.getStringValue(COLUMN_RESPONSE,cls));
			
			dataframe.setValue(COLUMN_FID,patient_id,patient_id);
			dataframe.setValue(COLUMN_IID,patient_id,patient_id);
			dataframe.setValue(COLUMN_CLASS,patient_id,cls);
			dataframe.setValue(COLUMN_SEX,patient_id,sex);
			dataframe.setValue(COLUMN_AGE,patient_id,age);
			dataframe.setValue(COLUMN_TREATMENT,patient_id,treatment);
			dataframe.setValue(COLUMN_NUMTREATMENTS,patient_id,numtreatments);
			dataframe.setValue(COLUMN_RESPONSE,patient_id,response);
		}
		CFileHelper.writeFile(outfile,dataframe.toString());
	}

	public void createBfile()
	{
		//plink --tfile gwas --make-bed
		System.out.println("creating bfile from tfile and tfam entries");
		CCommandLine commands=new CCommandLine("plink");
		commands.setWorkingDir(params.outdir);
		commands.addArg("--tfile",params.name);
		commands.addArg("--make-bed");
		commands.addArg("--out",params.name);
		IRuntimeService runtimeService=new RuntimeServiceImpl();
		runtimeService.exec(commands);
	}
	
	public static class Params
	{
		public String dir="c:\\projects\\analysis\\GWAS\\GWAS_FinalReport\\";
		public String snpfile="Human610-Quadv1_B.txt";
		public String genofile="FinalReport_extr.csv";
		public String phenofile="hcv.phe";
		public String snplistfile=null;
		//public String samplefile="sample_sheet.txt";
		//public String classfile="classes.txt";
		public String outdir="c:\\temp\\";
		public String name="gwas";
		public boolean rsonly=true;
		public boolean minimize=true;
		
		public Params(){}
		
		public Params(String[] args)
		{
			setParams(args);
		}
		
		public void setParams(String[] argv)
		{
			Options options = new Options();
			options.addOption("dir", true,"input directory (c:\\projects\\analysis\\GWAS\\GWAS_FinalReport\\)");
			options.addOption("snpfile", true,"contains information about each SNP (Human610-Quadv1_B.txt)");
			options.addOption("genofile", true,"contains genotyping data for each sample (FinalReport_extr.csv)");
			options.addOption("phenofile", true,"phenotype file (hcv.phe)");
			options.addOption("snplistfile", true," list of snps to include (hcv.snplist)");
			options.addOption("outdir", true,"output directory (c:\\temp\\)");
			options.addOption("name", true,"output directory (gwas)");
			options.addOption("rsonly", true,"Include only reference SNPs (rs prefix) (true)");
			options.addOption("minimize", true,"Include genotypes for patients in pheno file only");
			
			CommandLine line=PlinkHelper.parseCommandLine(argv,options);
			if (line.hasOption("dir")) this.dir=CFileHelper.normalizeDirectory(line.getOptionValue("dir"));
			if (line.hasOption("snpfile")) this.snpfile=fixPath(line.getOptionValue("snpfile"));
			if (line.hasOption("genofile")) this.genofile=fixPath(line.getOptionValue("genofile"));
			if (line.hasOption("phenofile")) this.phenofile=fixPath(line.getOptionValue("phenofile"));
			if (line.hasOption("snplistfile")) this.snplistfile=fixPath(line.getOptionValue("snplistfile"));
			if (line.hasOption("outdir")) this.outdir=CFileHelper.normalizeDirectory(line.getOptionValue("outdir"));
			if (line.hasOption("name")) this.name=line.getOptionValue("name");
			if (line.hasOption("rsonly")) this.rsonly=Boolean.valueOf(line.getOptionValue("rsonly"));
			if (line.hasOption("minimize")) this.minimize=Boolean.valueOf(line.getOptionValue("minimize"));
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
