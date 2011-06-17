package org.vardb.util.tools.snps;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.vardb.util.CDataFrame;
import org.vardb.util.CException;
import org.vardb.util.CFileHelper;
import org.vardb.util.CStringHelper;
import org.vardb.util.CTable;

public final class SnpReportGenerator
{
	private SnpReportGenerator(){}
	
	private static final String ASSOC_SUFFIX=".assoc";
	private static final String ASSOC_ADJUSTED_SUFFIX=".assoc.adjusted";
	
	private static final String ASSOC_LINEAR_SUFFIX=".assoc.linear";
	private static final String ASSOC_LINEAR_ADJUSTED_SUFFIX=".assoc.linear.adjusted";
	
	private static final String QASSOC_SUFFIX=".qassoc";
	private static final String QASSOC_ADJUSTED_SUFFIX=".qassoc.adjusted";
	private static final String QASSOC_MEANS_SUFFIX=".qassoc.means";
	
	public static void createReport(Params params)
	{
		System.out.println("extracting SNPs from directory "+params.indir);
		CFileHelper.createDirectory(params.reportdir);
		Map<String,Collection<SnpData>> sharedSnps=new HashMap<String,Collection<SnpData>>();
		addSharedSnps(sharedSnps,extractSnps(params,AnalysisType.ASSOC));
		addSharedSnps(sharedSnps,extractSnps(params,AnalysisType.QASSOC));
		addSharedSnps(sharedSnps,extractSnps(params,AnalysisType.ASSOC_LINEAR));	
		addFrequencyData(params,sharedSnps);
		reportSharedSnps(params,sharedSnps);
		createAnalysisDataFrames(params,sharedSnps);
	}

	private static void createAnalysisDataFrames(Params params, Map<String,Collection<SnpData>> sharedSnps)
	{
		for (String identifier : sharedSnps.keySet())
		{
			createAnalysisDataFrame(params, identifier, sharedSnps.get(identifier));
		}
	}
	
	private static void createAnalysisDataFrame(Params params, String identifier, Collection<SnpData> snps)
	{

		System.out.println("trying to read phenotype file "+params.plinkdir+params.phenofile);
		CDataFrame dataframe=CDataFrame.parseTabFile(params.plinkdir+params.phenofile);
		Set<String> ids=getSnpIds(snps);
		if (!ids.isEmpty())
		{
			GenotypeExtractor extractor=new GenotypeExtractor(new GenotypeExtractor.Params(params.plinkdir,params.bfile));
			CDataFrame genotypes=extractor.extract(ids,params.indicators);
			dataframe.appendColumns(genotypes);
		}
		else System.out.println("no SNPs found for "+identifier);		
		String outfile=params.reportdir+identifier+".data";
		System.out.println("writing phenotype+SNP file: "+outfile);
		PlinkHelper.writeDataFile(outfile,dataframe);
	}
	
	private static void addFrequencyData(Params params, Map<String,Collection<SnpData>> sharedSnps)
	{
		Set<String> ids=getSharedSnpIds(sharedSnps);
		Map<String,PlinkHelper.SnpFrequencyData> freqs=PlinkHelper.getSnpFrequencies(params.plinkdir+params.freqfile,ids);
		for (Collection<SnpData> snps : sharedSnps.values())
		{
			for (SnpData snp : snps)
			{
				PlinkHelper.SnpFrequencyData freq=freqs.get(snp.snp);
				snp.setFrequencyData(freq);
			}
		}
	}
	
	private static Set<String> getSharedSnpIds(Map<String,Collection<SnpData>> sharedSnps)
	{
		Set<String> ids=new HashSet<String>();
		for (String identifier : sharedSnps.keySet())
		{
			Collection<SnpData> snps=sharedSnps.get(identifier);
			for (SnpData snp : snps)
			{
				if (!ids.contains(snp.snp))
					ids.add(snp.snp);
			}
		}
		return ids;
	}
	
	private static Set<String> getSnpIds(Collection<SnpData> snps)
	{
		Set<String> ids=new HashSet<String>();
		for (SnpData snp : snps)
		{
			if (!ids.contains(snp.snp))
				ids.add(snp.snp);
		}
		return ids;
	}
	
	///////////////////////////////////////////////////////
	
	public static Map<String,Collection<SnpData>> extractSnps(Params params, AnalysisType type)
	{
		System.out.println("extracting SNPs of type "+type);
		Map<String,Collection<SnpData>> map=new HashMap<String,Collection<SnpData>>();
		List<String> identifiers=getIdentifiers(params.indir,type);
		for (String identifier : identifiers)
		{
			Collection<SnpData> topsnps=getSnpData(params,identifier,type);
			writeSnps(params,type,topsnps,identifier);
			map.put(identifier,topsnps);
		}
		return map;
	}
	
	private static Collection<SnpData> getSnpData(Params params, String identifier, AnalysisType type)
	{
		Map<String,SnpData> topsnps=getTopSnps(params,identifier,type);
		addSnpInfo(params,identifier,type,topsnps);
		addSnpMeans(params,identifier,type,topsnps);
		return topsnps.values();
	}

	private static Map<String,SnpData> getTopSnps(Params params, String identifier, AnalysisType type)
	{
		Map<String,SnpData> topsnps=new LinkedHashMap<String,SnpData>();
		String filename=params.indir+identifier+type.getAdjustedSuffix();
		Scanner scanner=CFileHelper.createScanner(filename,"\n");
		scanner.nextLine(); //skip the header line
		while(scanner.hasNext())
		{
			String line=scanner.nextLine();
			String[] tabs=line.split("\\s+");
			SnpData snp=type.createInstance();
			snp.setAdjusted(tabs);
			if (snp.getBonf()>=params.cutoff)
				break;
			topsnps.put(snp.snp,snp);
		}
		return topsnps;
	}
	
	private static void addSnpInfo(Params params, String identifier, AnalysisType type, Map<String,SnpData> topsnps)
	{
		String filename=params.indir+identifier+type.getSuffix();
		Scanner scanner=CFileHelper.createScanner(filename,"\n");
		scanner.nextLine(); //skip the header line
		while(scanner.hasNext())
		{
			String line=scanner.nextLine();
			String[] tabs=line.split("\\s+");
			String snpID=tabs[2];
			if (!topsnps.containsKey(snpID))
				continue;
			SnpData snp=topsnps.get(snpID);
			snp.setInfo(tabs);
		}
	}
	
	private static void addSnpMeans(Params params, String identifier, AnalysisType analysisType, Map<String,SnpData> topsnps)
	{
		if (analysisType.getMeansSuffix()==null)
			return;
		String filename=params.indir+identifier+analysisType.getMeansSuffix();
		Scanner scanner=CFileHelper.createScanner(filename,"\n");
		scanner.nextLine(); //skip the header line
		while(scanner.hasNext())
		{
			String line=scanner.nextLine();
			String[] tabs=line.split("\\s+");
			String snpID=tabs[2];
			if (!topsnps.containsKey(snpID))
				continue;
			QAssocSnpData snp=(QAssocSnpData)topsnps.get(snpID);
			snp.setMeans(tabs);
		}
	}
	
	private static List<String> getIdentifiers(String dir, AnalysisType type)
	{
		List<String> filenames=CFileHelper.listFiles(dir,type.getSuffix());
		List<String> identifiers=new ArrayList<String>();
		for (String filename : filenames)
		{
			String name=CFileHelper.getIdentifierFromFilename(filename, type.getSuffix());
			identifiers.add(name);			
		}
		return identifiers;
	}
	
	private static void writeSnps(Params params, AnalysisType type, Collection<SnpData> snps, String identifier)
	{
		CTable table=type.createInstance().createTable(snps);
		String tracks=type.createInstance().createTracks(identifier,snps);
		CFileHelper.writeFile(params.reportdir+identifier+"."+type.name().toLowerCase()+".report",table.toString());
		CFileHelper.writeFile(params.reportdir+identifier+"."+type.name().toLowerCase()+".tracks",tracks);
	}
	
	private static String findMainChromosome(Collection<? extends ISnpData> snps)
	{
		Map<String,Integer> counts=new HashMap<String,Integer>();
		for (ISnpData item : snps)
		{
			String chr=item.getChr();
			Integer count=counts.get(chr);
			if (count==null)
				count=1;
			else count++;
			counts.put(chr,count);
		}
		int max=0;
		String maxchr=null;
		for (String chr : counts.keySet())
		{
			int count=counts.get(chr);
			if (count>max)
			{
				max=count;
				maxchr=chr;
			}
		}
		return maxchr;
	}
	
	public static class AssocSnpData extends SnpData
	{
		public String a1;
		public String f_a;
		public String f_u;
		public String a2;
		public String chisq;
		public String p;
		public String or;

		//CHR         SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR
		//0 MitoT10239C      10239    B 0.004902     0.05    A        13.51    0.0002367       0.0936
		public void setInfo(String[] tabs)
		{
			bp=Integer.valueOf(tabs[3]);
			a1=tabs[4];
			f_a=tabs[5];
			f_u=tabs[5];
			a2=tabs[5];
			chisq=tabs[5];
			p=tabs[5];
			or=tabs[5];
		}
		
		@Override
		protected void addInfoColumns(CTable table)
		{
			table.getHeader().add("a1");
			table.getHeader().add("f_a");
			table.getHeader().add("f_u");
			table.getHeader().add("a2");
			table.getHeader().add("chisq");
			table.getHeader().add("p");
			table.getHeader().add("or");
		}
		
		@Override
		protected void addInfoColumns(CTable.Row row)
		{
			row.add("a1");
			row.add("f_a");
			row.add("f_u");
			row.add("a2");
			row.add("chisq");
			row.add("p");
			row.add("or");
		}
	}
	
	public static class AssocLinearSnpData extends SnpData
	{
		public String a1;
		public String test;
		public Integer nmiss;
		public String beta;
		public String stat;
		public String p;
		
		//CHR         SNP         BP   A1       TEST    NMISS       BETA         STAT            P
	   	//0 MitoT10239C      10239    B        ADD      237     0.6408       0.5264       0.5991
		public void setInfo(String[] tabs)
		{
			bp=Integer.valueOf(tabs[3]);
			a1=tabs[4];
			test=tabs[5];
			nmiss=Integer.parseInt(tabs[6]);
			beta=tabs[7];
			stat=tabs[8];
			p=tabs[9];
		}
		
		@Override
		protected void addInfoColumns(CTable table)
		{
			table.getHeader().add("a1");
			table.getHeader().add("test");
			table.getHeader().add("nmiss");
			table.getHeader().add("beta");
			table.getHeader().add("stat");
			table.getHeader().add("p");
		}
		
		@Override
		protected void addInfoColumns(CTable.Row row)
		{
			row.add(a1);
			row.add(test);
			row.add(nmiss);
			row.add(beta);
			row.add(stat);
			row.add(p);
		}
	}
	
	public static class QAssocSnpData extends SnpData
	{
		public Integer nmiss;
		public String beta;
		public String se;
		public String r2;
		public String t;
		public String p;		
		
		public SnpMean g11=new SnpMean();
		public SnpMean g12=new SnpMean();
		public SnpMean g22=new SnpMean();
		
		//CHR         SNP         BP    NMISS       BETA         SE         R2        T            P
		//0 MitoT10239C      10239      284    -0.3301     0.1449    0.01807   -2.278      0.02347
		public void setInfo(String[] tabs)
		{
			bp=Integer.valueOf(tabs[3]);
			nmiss=Integer.valueOf(tabs[4]);
			beta=tabs[5];
			se=tabs[6];
			r2=tabs[7];
			t=tabs[8];
			p=tabs[9];
		}
		
		//CHR         SNP  VALUE      G11      G12      G22
		//0 MitoT10239C   GENO      B/B      B/A      A/A
		public void setMeans(String[] tabs)
		{
			MeanType type=MeanType.valueOf(tabs[3]);
			g11.setValue(type,tabs[4]);
			g12.setValue(type,tabs[5]);
			g22.setValue(type,tabs[6]);
		}
		
		@Override
		protected void addInfoColumns(CTable table)
		{
			table.getHeader().add("nmiss");
			table.getHeader().add("beta");
			table.getHeader().add("se");
			table.getHeader().add("r2");
			table.getHeader().add("t");
			table.getHeader().add("p");
		}
		
		protected void addInfoColumns(CTable.Row row)
		{
			row.add(nmiss);
			row.add(beta);
			row.add(se);
			row.add(r2);
			row.add(t);
			row.add(p);
		}
		
		@Override
		protected void addMeanColumns(CTable table)
		{
			table.getHeader().add("G11counts");
			table.getHeader().add("G12counts");
			table.getHeader().add("G22counts");
			
			table.getHeader().add("G11freq");
			table.getHeader().add("G12freq");
			table.getHeader().add("G22freq");
			
			table.getHeader().add("G11mean");
			table.getHeader().add("G12mean");
			table.getHeader().add("G22mean");
			
			table.getHeader().add("G11sd");
			table.getHeader().add("G12sd");
			table.getHeader().add("G22sd");
		}

		@Override
		protected void addMeanColumns(CTable.Row row)
		{
			row.add(g11.counts);
			row.add(g12.counts);
			row.add(g22.counts);
			
			row.add(g11.freq);
			row.add(g12.freq);
			row.add(g22.freq);
			
			row.add(g11.mean);
			row.add(g12.mean);
			row.add(g22.mean);
			
			row.add(g11.sd);
			row.add(g12.sd);
			row.add(g22.sd);
		}
	}

	public abstract static class SnpData implements ISnpData
	{
		public String chr;
		public String snp;
		public Integer bp;
		
		public String a1;
		public String a2;
		public String maf;
		
		public String unadj;	
		public String gc;
		public String qq;
		public String bonf;
		public String holm;
		public String sidakSs;
		public String sidakSd;
		public String fdrBh;
		public String fdrBy;

		//CHR         SNP      UNADJ         GC         QQ       BONF       HOLM   SIDAK_SS   SIDAK_SD     FDR_BH     FDR_BY
		//CHR         SNP      UNADJ         GC       BONF       HOLM   SIDAK_SS   SIDAK_SD     FDR_BH     FDR_BY
	  	//1   rs2774941 2.838e-006 2.838e-006          1          1     0.7357     0.7357     0.1901          1
		public void setAdjusted(String[] tabs)
		{
			chr=tabs[1];
			snp=tabs[2];
			unadj=tabs[3];
			gc=tabs[4];
			qq=tabs[5];
			bonf=tabs[6];
			holm=tabs[7];
			sidakSs=tabs[8];
			sidakSd=tabs[9];
			fdrBh=tabs[10];
			fdrBy=tabs[11];
		}
		
		public abstract void setInfo(String[] tabs);
		
		public double getBonf()
		{
			if (bonf.equals("INF"))
				return 1;
			return Double.parseDouble(bonf);
		}
		
		protected void addInfoColumns(CTable table){}

		protected void addMeanColumns(CTable table){}
		
		protected void addInfoColumns(CTable.Row row){}
		
		protected void addMeanColumns(CTable.Row row){}
		
		public CTable createTable(Collection<SnpData> snps)
		{
			CTable table=new CTable();
			table.getHeader().add("chr");
			table.getHeader().add("snp");
			table.getHeader().add("bp");			
			addInfoColumns(table);			
			table.getHeader().add("gc");
			table.getHeader().add("bonf");
			table.getHeader().add("holm");
			table.getHeader().add("sidakSs");
			table.getHeader().add("sidakSd");
			table.getHeader().add("fdrBh");
			table.getHeader().add("fdrBy");
			table.getHeader().add("a1");
			table.getHeader().add("a2");
			table.getHeader().add("maf");
			addMeanColumns(table);			
			for (SnpData item : snps)
			{
				item.addRow(table);
			}			
			return table;
		}

		public void addRow(CTable table)
		{
			CTable.Row row=table.addRow();
			row.add(chr);
			row.add(snp);
			row.add(bp);		
			addInfoColumns(row);			
			row.add(gc);
			row.add(bonf);
			row.add(holm);
			row.add(sidakSs);
			row.add(sidakSd);
			row.add(fdrBh);
			row.add(fdrBy);
			row.add(a1);
			row.add(a2);
			row.add(maf);
			addMeanColumns(row);
		}
		
		private String createTracks(String name, Collection<SnpData> snps)
		{
			StringBuilder buffer=new StringBuilder();			
			buffer.append("browser position chr").append(findMainChromosome(snps)).append("\n");
			buffer.append("track name="+name+" description=\""+name+"\" visibility=2\n");			
			for (SnpData item : snps)
			{
				item.addTrack(buffer);
			}
			return buffer.toString();
		}	
		
		private void addTrack(StringBuilder buffer)
		{
			String chr=this.chr.equals("23") ? "X" : this.chr;
			buffer.append("chr").append(chr).append(" ");
			buffer.append(bp).append(" ");
			buffer.append(bp).append(" ");
			buffer.append(snp).append("\n");
		}
		
		public String getChr()
		{
			if ("23".equals(chr))
				return "X";
			return chr;
		}
		
		public void setFrequencyData(PlinkHelper.SnpFrequencyData freq)
		{
			if (freq==null)
				throw new CException("Frequency data is null for SNP "+snp);
			if (!this.snp.equals(freq.snp))
				throw new CException("SNP names do not match: this.snp="+snp+", freq.snp="+freq.snp);
			this.a1=freq.a1;
			this.a2=freq.a2;
			this.maf=freq.maf;
		}
	}
	
	public enum MeanType{GENO,COUNTS,FREQ,MEAN,SD};
	
	public static class SnpMean
	{
		public String geno;
		public String counts;
		public String freq;
		public String mean;
		public String sd;
		
		public void setValue(MeanType type, String value)
		{
			switch(type)
			{
			case GENO:
				geno=value;
				break;
			case COUNTS:
				counts=value;
				break;
			case FREQ:
				freq=value;
				break;
			case MEAN:
				mean=value;
				break;
			case SD:
				sd=value;
				break;
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////////
	
	private static void addSharedSnps(Map<String,Collection<SnpData>> sharedSnps, Map<String,Collection<SnpData>> snps)
	{
		for (String identifier : snps.keySet())
		{
			String root=findRoot(identifier);
			if (!sharedSnps.containsKey(root))
				sharedSnps.put(root,new ArrayList<SnpData>());
			Collection<SnpData> collection=sharedSnps.get(root);
			collection.addAll(snps.get(identifier));
		}
	}
	
	private static String findRoot(String identifier)
	{
		int index=identifier.indexOf('_');
		if (index==-1)
		{
			System.out.println("NOTE: no root in identifier name: "+identifier);
			return identifier;//throw new CException("can't find root in identifier "+identifier);
		}
		return identifier.substring(0,index);
	}
	
	private static void reportSharedSnps(Params params, Map<String,Collection<SnpData>> sharedSnps)
	{
		System.out.println("generating reports for shared snps");
		for (String identifier : sharedSnps.keySet())
		{
			reportSharedSnps(params,identifier,sharedSnps.get(identifier));
		}
	}
	
	private static void reportSharedSnps(Params params, String identifier, Collection<SnpData> snps)
	{
		System.out.println("generating shared snp report for "+identifier);
		Map<String,SnpCountData> countmap=new HashMap<String,SnpCountData>();
		for (SnpData snpdata : snps)
		{
			if (!countmap.containsKey(snpdata.snp))
				countmap.put(snpdata.snp,new SnpCountData(snpdata));
			SnpCountData count=countmap.get(snpdata.snp); 
			count.add(snpdata);
		}
	
		List<SnpCountData> counts=SnpCountData.sortCountData(countmap.values());
		CTable report=SnpCountData.getReport(counts);
		CFileHelper.writeFile(params.reportdir+identifier+".report",report.toString());
		String tracks=SnpCountData.getTracks(identifier,counts);
		CFileHelper.writeFile(params.reportdir+identifier+".tracks",tracks);
	}
	
	private static class SnpCountData implements ISnpData
	{
		public String snp;
		public String chr;
		public Integer bp;
		public Integer count=0;
		public String a1;
		public String a2;
		public String maf;
		public String nchrobs;
		public List<Double> pvalues=new ArrayList<Double>();
		
		public SnpCountData(SnpData snpdata)
		{
			this.snp=snpdata.snp;
			this.chr=snpdata.chr;
			this.bp=snpdata.bp;
			this.a1=snpdata.a1;
			this.a2=snpdata.a2;
			this.maf=snpdata.maf;
		}
		
		public void add(SnpData data)
		{
			count=count+1;
			pvalues.add(data.getBonf());
		}
		
		public void sort()
		{
			Collections.sort(pvalues);
		}
		
		public static CTable getReport(Collection<SnpCountData> counts)
		{
			CTable table=new CTable();
			table.getHeader().add("snp");
			table.getHeader().add("chr");
			table.getHeader().add("bp");
			table.getHeader().add("count");
			table.getHeader().add("pvalues");
			table.getHeader().add("a1");
			table.getHeader().add("a2");
			table.getHeader().add("maf");
			table.getHeader().add("nchrobs");
			for (SnpCountData count : counts)
			{
				count.getRow(table);
			}
			return table;
		}
		
		private void getRow(CTable table)
		{
			CTable.Row row=table.addRow();
			row.add(snp);
			row.add(chr);
			row.add(bp);
			row.add(count);
			row.add(CStringHelper.join(pvalues,", "));
			row.add(a1);
			row.add(a2);
			row.add(maf);
			row.add(nchrobs);
		}
		
		public static List<SnpCountData> sortCountData(Collection<SnpCountData> coll)
		{
			List<SnpCountData> counts=new ArrayList<SnpCountData>();
			counts.addAll(coll);
			for (SnpCountData count : counts)
			{
				count.sort();
			}		
			Collections.sort(counts,new SnpCountDataComparator());
			return counts;
		}
		
		public static String getTracks(String identifier, Collection<SnpCountData> coll)
		{
			StringBuilder buffer=new StringBuilder();			
			buffer.append("browser position chr").append(findMainChromosome(coll)).append("\n");
			buffer.append("track name="+identifier+" description=\""+identifier+"\" visibility=2\n");
			for (SnpCountData countdata : coll)
			{
				countdata.addTrack(buffer);
			}
			return buffer.toString();
		}
		
		private void addTrack(StringBuilder buffer)
		{
			String chr=this.chr.equals("23") ? "X" : this.chr;
			buffer.append("chr").append(chr).append(" ");
			buffer.append(bp).append(" ");
			buffer.append(bp).append(" ");
			buffer.append(snp).append("\n");
		}
				
		public String getChr()
		{
			if ("23".equals(chr))
				return "X";
			return chr;
		}
	}

	private static class SnpCountDataComparator implements Comparator<SnpCountData>
	{
		@Override
		public int compare(SnpCountData obj1, SnpCountData obj2)
		{
			if (obj2.count.equals(obj1.count)) // if the counts are the same sort by p-value
				return obj1.pvalues.get(0).compareTo(obj2.pvalues.get(0));
			return obj2.count.compareTo(obj1.count);
		}
	}
	
	public interface ISnpData
	{
		public String getChr();
	}
	
	public enum AnalysisType
	{
		ASSOC(QAssocSnpData.class,ASSOC_SUFFIX,ASSOC_ADJUSTED_SUFFIX),
		ASSOC_LINEAR(AssocLinearSnpData.class,ASSOC_LINEAR_SUFFIX,ASSOC_LINEAR_ADJUSTED_SUFFIX),
		QASSOC(QAssocSnpData.class,QASSOC_SUFFIX,QASSOC_ADJUSTED_SUFFIX,QASSOC_MEANS_SUFFIX);
		
		private Class<? extends SnpData> cls;
		private String suffix;
		private String adjusted;
		private String means;
		
		AnalysisType(Class<? extends SnpData> cls, String suffix, String adjusted)
		{
			this.cls=cls;
			this.suffix=suffix;
			this.adjusted=adjusted;

		}
		
		AnalysisType(Class<? extends SnpData> cls, String suffix, String adjusted, String means)
		{
			this(cls,suffix,adjusted);
			this.means=means;
		}
		
		public String getSuffix(){return suffix;}
		public String getAdjustedSuffix(){return adjusted;}
		public String getMeansSuffix(){return means;}
		
		public SnpData createInstance()
		{
			try
			{
				return cls.newInstance();
			}
			catch(Exception e)
			{
				throw new CException(e);
			}
		}
	}
	
	public static class Params
	{
		public String plinkdir;
		public String bfile;
		public String indir;
		public String outdir;
		public double cutoff=0.2;
		public String phenofile;
		public String freqfile;
		public String reportdir;
		public boolean indicators=true;
		
		public Params(String plinkdir, String bfile, String outdir)
		{
			this.plinkdir=plinkdir;
			this.bfile=bfile;
			this.outdir=outdir;
			init();
		}
		
		public Params(String[] args)
		{
			setParams(args);
		}
		
		public void setParams(String[] args)
		{
			Options options = new Options();
			options.addOption("plinkdir", true, "input directory (c:\\projects\\analysis\\GWAS\\plink\\)");
			options.addOption("bfile", true, "root name for bed/bim/fam files");
			options.addOption("indir", true, "output directory (c:\\projects\\analysis\\GWAS\\plink\\)");
			options.addOption("outdir", true, "output directory (c:\\projects\\analysis\\GWAS\\plink\\)");
			options.addOption("cutoff", true, "threshold for inclusion of SNPs in analysis");
			
			CommandLine line=PlinkHelper.parseCommandLine(args,options);
	        if (line.hasOption("plinkdir")) this.plinkdir=line.getOptionValue("plinkdir");
	        if (line.hasOption("bfile")) this.bfile=line.getOptionValue("bfile");
	        if (line.hasOption("indir")) this.indir=line.getOptionValue("indir");
	        if (line.hasOption("outdir")) this.outdir=line.getOptionValue("outdir");
	        if (line.hasOption("cutoff")) this.cutoff=Double.valueOf(line.getOptionValue("cutoff"));	        
	        init();	        
		}
		
		private void init()
		{
			this.plinkdir=CFileHelper.normalizeDirectory(this.plinkdir);
			this.indir=CFileHelper.normalizeDirectory(this.indir);
			this.outdir=CFileHelper.normalizeDirectory(this.outdir);
			this.phenofile=bfile+".phe";
			this.freqfile=bfile+".frq";
			this.reportdir=outdir+"reports/";
		}
		
		@Override
		public String toString()
		{
			return CStringHelper.toString(this);
		}
	}
	
	public static void main(String[] args)
	{
		//Params params=new Params(args[0],args[1],args[2]);
		//if (args.length>3)
		//	params.cutoff=Double.parseDouble(args[3]);

		Params params=new Params(args);
		System.out.println("params: "+params);
		
		//PlinkHelper.getSnpFrequencies(params.plinkdir+params.bfile+".frq",Arrays.asList("rs8099917","rs12979860"));
		SnpReportGenerator.createReport(params);
	}
}
