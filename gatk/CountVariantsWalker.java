import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.IntervalBinding;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

/**
 * Prints the base counts at each position
 * mvn install && java -Xmx2g -cp ./target/vardb-util-0.1.0.BUILD-SNAPSHOT.jar;../lib/GenomeAnalysisTK.jar org.broadinstitute.sting.gatk.CommandLineGATK -T CountVariants -R U:\ref\HCV-KT9.fasta -I U:\OLD2bam\nextgen3-2G__HCV-KT9.bam -dt NONE -O U:\counts.txt
 * java -Xmx10g -cp ../target/gatk-walkers.jar;../../lib/GenomeAnalysisTK.jar org.broadinstitute.sting.gatk.CommandLineGATK -T CountVariants -R U:\ref\HCV-KT9.fasta -I U:\bam\nextgen1-1A__HCV-KT9.bam -dt NONE --ntcounts U:\ntcounts.txt --codoncounts U:\codoncounts.txt --aacounts U:\aacounts.txt -L HCV-KT9:3420-5312 -L HCV-KT9:6258-7598
 * java -Xmx10g -cp ../target/gatk-walkers.jar;../../lib/GenomeAnalysisTK.jar org.broadinstitute.sting.gatk.CommandLineGATK -T CountVariants -R U:\ref\HCV-KT9.fasta -I U:\bam\CTE247-21__HCV-KT9.bam -dt NONE --ntcounts U:\ntcounts.txt --codoncounts U:\codoncounts.txt --aacounts U:\aacounts.txt -L HCV-KT9:3420-5312 -L HCV-KT9:6258-7598
 * 
 * 
 */
//@By(DataSource.REFERENCE)
public class CountVariantsWalker extends LocusWalker<Integer, Integer> implements TreeReducible<Integer>
{
	//@Input(fullName = "input_file", shortName = "I", doc = "SAM or BAM file(s)", required = false)
	
	@Output(fullName="ntcounts", shortName="o1", required=true)
	protected PrintStream out1;
	
	@Output(fullName="codoncounts", shortName="o2", required=true)
	protected PrintStream out2;
	
	@Output(fullName="aacounts", shortName="o3", required=true)
	protected PrintStream out3;
	
	protected Integer baseQuality=30;
	protected Integer mapQuality=30;	
	protected List<Region> regions=new ArrayList<Region>();
	
	static String TAB="\t";
	
	public void initialize()
	{
		initRegions();
		writeBaseCountHeader(out1);
		writeCodonCountsHeader(out2);
		writeAminoAcidCountsHeader(out3);
	}

	public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context)
    {
    	System.out.println("position="+context.getPosition());
    	Region region=getRegion(context.getPosition());
    	ReadBackedPileup pileup = context.getBasePileup().getBaseAndMappingFilteredPileup(baseQuality,mapQuality);
    	if (region.addBases(context,pileup))
    	{
    		List<CodonCount> codoncounts=region.countCodons();
    		writeCodonCounts(out2,region,codoncounts);
    		List<AminoAcidCount> aacounts=region.countAminoAcids(codoncounts);
    		writeAminoAcidCounts(out3,region,aacounts);
    	}
    	writeBaseCounts(out1,ref,context,pileup);
    	return 1;
    }
	
	/////////////////////////////////////////

	private void writeCodonCountsHeader(PrintStream out)
	{
		out.println("position"+TAB+"aanum"+TAB+"aa"+TAB+"codon"+TAB+"rank"+TAB+"count"+TAB+"freq");
	}
	
	private void writeCodonCounts(PrintStream out, Region region, List<CodonCount> counts)
	{
		for (CodonCount count : counts)
		{
			//System.out.println("codon="+count.codon+"="+count.count);
			out.print(region.codonstart);
			out.print(TAB+region.aanum);
			out.print(TAB+count.codon.getAminoAcid().getCode());
			out.print(TAB+count.codon.name());
			out.print(TAB+count.rank);
			out.print(TAB+count.count);
			out.printf(TAB+"%2.3f",count.freq);
			out.println();
		}
		region.clear();
	}
	
	/////////////////////////////////////////
	
	private void writeAminoAcidCountsHeader(PrintStream out)
	{
		out.println("position"+TAB+"aanum"+TAB+"aa"+TAB+"rank"+TAB+"count"+TAB+"freq");
	}

	private void writeAminoAcidCounts(PrintStream out, Region region, List<AminoAcidCount> counts)
	{
		for (AminoAcidCount count : counts)
		{
			//System.out.println("aa="+count.aa+"="+count.count);
			out.print(region.codonstart);
			out.print(TAB+region.aanum);
			out.print(TAB+count.aa.getCode());
			out.print(TAB+count.rank);
			out.print(TAB+count.count);
			out.printf(TAB+"%2.3f",count.freq);
			out.println();
		}
		region.clear();
	}
	
	///////////////////////////////////////////
	
	private Region getRegion(long position)
	{
		for (Region region : regions)
		{
			if (region.inRegion(position))
				return region;
		}
		throw new RuntimeException("Region should not be null if using intervals: position="+position);
	}
	
	@Override
    public Integer reduceInit() { return 0; }
    
	@Override
    public Integer reduce(Integer value, Integer sum)
    {
        return treeReduce(sum,value);
    }
    
	@Override
    public Integer treeReduce(Integer lhs, Integer rhs)
    {
        return lhs + rhs;
    }

	private void initRegions()
	{
		for (IntervalBinding<Feature> interval : super.getToolkit().getArguments().intervals)
    	{
			//System.out.println("interval: "+interval.toString());
    		List<GenomeLoc> locations=interval.getIntervals(super.getToolkit());
    		for (GenomeLoc loc : locations)
    		{
    			regions.add(new Region(loc.getStart(),loc.getStop()));
    		}
    	}
	}

	private void writeBaseCountHeader(PrintStream out)
	{
		out.println("position"+TAB+"refnt"+TAB+"a"+TAB+"c"+TAB+"g"+TAB+"t"+TAB+"depth");
	}
	
	private void writeBaseCounts(PrintStream out, ReferenceContext ref, AlignmentContext context, ReadBackedPileup pileup)
	{
		int[] counts=pileup.getBaseCounts();
    	out1.printf("%d\t",context.getPosition());
    	out.printf("%s\t",(char)ref.getBase());
    	out.printf("%d\t%d\t%d\t%d\t",counts[0],counts[1],counts[2],counts[3]);//A, C, G, T
    	out.printf("%d%n",pileup.depthOfCoverage()); 
	}
}

class Region
{
	protected int start;
	protected int end;
	protected Map<Integer,Map<String,Character>> map=new HashMap<Integer,Map<String,Character>>();
	protected int codonstart;
	protected int aanum=0;
	
	public Region(int start, int end)
	{
		//System.out.println(" start="+start+" stop="+end);
		this.start=start;
		this.end=end;
		for (int offset=0;offset<=2;offset++)
		{
			map.put(offset,new HashMap<String,Character>());
		}
	}
	
	public void clear()
	{
		//System.out.println("clearing region");
		for (int offset=0;offset<=2;offset++)
		{
			map.get(offset).clear();
		}	
	}
	
	public boolean inRegion(long position)
	{
		return (position>=start && position<=end);				
	}
	
	public boolean isStart(long position)
	{
		return (position==start);
	}
	
	public boolean addBases(AlignmentContext context, ReadBackedPileup pileup)
	{
		int position=(int)context.getPosition();
		int relpos=position-start;
		int offset=relpos%3;
		if (offset==0)
		{
			codonstart=position;
			aanum++;
		}
		//System.out.println("position="+position+", relpos="+relpos+", offset="+offset);
		List<GATKSAMRecord> reads=pileup.getReads();
    	byte[] bases=pileup.getBases();
    	if (reads.size() != bases.length)
    		throw new RuntimeException("reads.size="+reads.size()+" is not equal to bases.length="+bases.length);
    	Map<String,Character> nts=map.get(offset);
    	for (int index=0; index<reads.size();index++)
    	{
    		GATKSAMRecord read=reads.get(index);
    		nts.put(read.getReadName(),(char)bases[index]);
    	}
    	return (offset==2);
	}
	
	public List<CodonCount> countCodons()
	{	
		Map<String,Character> nt1=map.get(0);
		Map<String,Character> nt2=map.get(1);
		Map<String,Character> nt3=map.get(2);
		Map<String,CodonCount> codons=new HashMap<String,CodonCount>();
		for (String name : nt1.keySet())
		{
			if (nt2.containsKey(name) && nt3.containsKey(name))
			{
				//System.out.println("name="+name+", base1="+base1+", base2="+base2+", base3="+base3);				
				String codon=""+nt1.get(name)+nt2.get(name)+nt3.get(name);
				CodonCount count=codons.get(codon);
				if (count==null)
					codons.put(codon,new CodonCount(codon));
				else count.increment();
			}
		}
		//System.out.println("*********************************************");
		List<CodonCount> counts=new ArrayList<CodonCount>();
		counts.addAll(codons.values());
		Collections.sort(counts,new CodonCountComparator());
		int rank=1;
		int sum=0;
		for (CodonCount count : counts)
		{
			count.rank=rank;
			rank++;
			sum+=count.count;
		}
		for (CodonCount count : counts)
		{
			count.freq=(float)count.count/(float)sum;
		}
		return counts;
	}
	
	public List<AminoAcidCount> countAminoAcids(List<CodonCount> codons)
	{	
		Map<AminoAcid,AminoAcidCount> aas=new HashMap<AminoAcid,AminoAcidCount>();
		for (CodonCount codon : codons)
		{
			AminoAcid aa=codon.codon.getAminoAcid();
			AminoAcidCount count=aas.get(aa);
			if (count==null)
				count=aas.put(aa,new AminoAcidCount(aa,codon.count));
			else count.increment(codon.count);
		}
		List<AminoAcidCount> counts=new ArrayList<AminoAcidCount>();
		counts.addAll(aas.values());
		Collections.sort(counts,new AminoAcidCountComparator());
		int rank=1;
		int sum=0;
		for (AminoAcidCount count : counts)
		{
			count.rank=rank;
			rank++;
			sum+=count.count;
		}
		for (AminoAcidCount count : counts)
		{
			count.freq=(float)count.count/(float)sum;
		}
		return counts;
	}
}

class CodonCount
{
	protected Codon codon;
	protected Integer count=1;
	protected Integer rank=null;
	protected Float freq=null;
	
	public CodonCount(String codon)
	{
		if (codon.length()!=3)
			throw new RuntimeException("codon must have length of 3: ["+codon+"]");
		this.codon=Codon.find(codon);
	}
	
	public void increment()
	{
		count++;
	}
}

class CodonCountComparator implements Comparator<CodonCount>
{
	public int compare(CodonCount a, CodonCount b)
	{
		return b.count.compareTo(a.count);
	}
}

class AminoAcidCount
{
	protected AminoAcid aa;
	protected Integer count=0;
	protected Integer rank=null;
	protected Float freq=null;
	
	public AminoAcidCount(AminoAcid aa, int codoncount)
	{
		this.aa=aa;
		increment(codoncount);
	}
	
	public void increment(int codoncount)
	{
		count+=codoncount;
	}
}

class AminoAcidCountComparator implements Comparator<AminoAcidCount>
{
	public int compare(AminoAcidCount a, AminoAcidCount b)
	{
		return b.count.compareTo(a.count);
	}
}

enum AminoAcid
{
	GLYCINE("G","Gly","Glycine"),
	ALANINE("A","Ala","Alanine"),
	SERINE("S","Ser","Serine"),
	THREONINE("T","Thr","Threonine"),
	CYSTEINE("C","Cys","Cysteine"),
	VALINE("V","Val","Valine"),
	LEUCINE("L","Leu","Leucine"),
	ISOLEUCINE("I","Ile","Isoleucine"),
	METHIONINE("M","Met","Methionine"),
	PROLINE("P","Pro","Proline"),
	PHENYLALANINE("F","Phe","Phenylalanine"),
	TYROSINE("Y","Tyr","Tyrosine"),
	TRYPTOPHAN("W","Trp","Tryptophan"),
	ASPARTIC_ACID("D","Asp","Aspartic Acid"),
	GLUTAMIC_ACID("E","Glu","Glutamic Acid"),
	ASPARAGINE("N","Asn","Asparagine"),
	GLUTAMINE("Q","Gln","Glutamine"),
	HISTIDINE("H","His","Histidine"),
	LYSINE("K","Lys","Lysine"),
	ARGININE("R","Arg","Arginine"),
	STOP("*","Stop","Stop");
	
	private final String code;
	private final String shortname;
	private final String longname;

	AminoAcid(String code, String shortname, String longname)
	{
		this.code=code;
		this.shortname=shortname;
		this.longname=longname;
	}
	
	public String getCode(){return this.code;}
	public String getShort(){return this.shortname;}
	public String getLongname(){return this.longname;}	

	public static AminoAcid find(char code)
	{
		return find(String.valueOf(code));
	}
	
	public static AminoAcid find(String code)
	{
		for (AminoAcid aa : AminoAcid.values())
		{
			if (aa.getCode().equals(code))
				return aa;
		}
		throw new RuntimeException("can't find AminoAcid ["+code+"]");
	}
}

enum Codon
{	
	ATT(AminoAcid.ISOLEUCINE),
	ATC(AminoAcid.ISOLEUCINE),
	ATA(AminoAcid.ISOLEUCINE),
	CTT(AminoAcid.LEUCINE),
	CTC(AminoAcid.LEUCINE),
	CTA(AminoAcid.LEUCINE),
	CTG(AminoAcid.LEUCINE),
	TTA(AminoAcid.LEUCINE),
	TTG(AminoAcid.LEUCINE),
	GTT(AminoAcid.VALINE),
	GTC(AminoAcid.VALINE),
	GTA(AminoAcid.VALINE),
	GTG(AminoAcid.VALINE),
	TTT(AminoAcid.PHENYLALANINE),
	TTC(AminoAcid.PHENYLALANINE),
	ATG(AminoAcid.METHIONINE),
	TGT(AminoAcid.CYSTEINE),
	TGC(AminoAcid.CYSTEINE),
	GCT(AminoAcid.ALANINE),
	GCC(AminoAcid.ALANINE),
	GCA(AminoAcid.ALANINE),
	GCG(AminoAcid.ALANINE),
	GGT(AminoAcid.GLYCINE),
	GGC(AminoAcid.GLYCINE),
	GGA(AminoAcid.GLYCINE),
	GGG(AminoAcid.GLYCINE),
	CCT(AminoAcid.PROLINE),
	CCC(AminoAcid.PROLINE),
	CCA(AminoAcid.PROLINE),
	CCG(AminoAcid.PROLINE),
	ACT(AminoAcid.THREONINE),
	ACC(AminoAcid.THREONINE),
	ACA(AminoAcid.THREONINE),
	ACG(AminoAcid.THREONINE),
	TCT(AminoAcid.SERINE),
	TCC(AminoAcid.SERINE),
	TCA(AminoAcid.SERINE),
	TCG(AminoAcid.SERINE),
	AGT(AminoAcid.SERINE),
	AGC(AminoAcid.SERINE),
	TAT(AminoAcid.TYROSINE),
	TAC(AminoAcid.TYROSINE),
	TGG(AminoAcid.TRYPTOPHAN),
	CAA(AminoAcid.GLUTAMINE),
	CAG(AminoAcid.GLUTAMINE),
	AAT(AminoAcid.ASPARAGINE),
	AAC(AminoAcid.ASPARAGINE),
	CAT(AminoAcid.HISTIDINE),
	CAC(AminoAcid.HISTIDINE),
	GAA(AminoAcid.GLUTAMIC_ACID),
	GAG(AminoAcid.GLUTAMIC_ACID),
	GAT(AminoAcid.ASPARTIC_ACID),
	GAC(AminoAcid.ASPARTIC_ACID),
	AAA(AminoAcid.LYSINE),
	AAG(AminoAcid.LYSINE),
	CGT(AminoAcid.ARGININE), 
	CGC(AminoAcid.ARGININE), 
	CGA(AminoAcid.ARGININE), 
	CGG(AminoAcid.ARGININE), 
	AGA(AminoAcid.ARGININE), 
	AGG(AminoAcid.ARGININE),
	TAA(AminoAcid.STOP),
	TAG(AminoAcid.STOP),
	TGA(AminoAcid.STOP);
	
	private final AminoAcid aminoAcid;

	Codon(AminoAcid aminoAcid)
	{
		this.aminoAcid=aminoAcid;
	}
	
	public AminoAcid getAminoAcid(){return this.aminoAcid;}
	
	public boolean isStopCodon(){return (this.aminoAcid==AminoAcid.STOP);}
	
	public String getRna()
	{
		return name().replace('T','U');
	}
	
	public static Codon find(String value)
	{
		value=value.toUpperCase();
		for (Codon codon : Codon.values())
		{
			if (codon.name().equals(value) || codon.getRna().equals(value))
				return codon;
		}
		return null;
	}
	
	public static List<Codon> getCodons(AminoAcid aa)
	{
		List<Codon> list=new ArrayList<Codon>();
		for (Codon codon : Codon.values())
		{
			if (codon.getAminoAcid()==aa)
				list.add(codon);
		}
		return list;
	}
	
	public String getNucleotide(int position)
	{
		if (position<0 || position>2)
			throw new RuntimeException("codon position should be between 0 and 2: "+position);
		return name().substring(position,position+1);
	}
	
	public static boolean isStopCodon(String value)
	{
		value=value.toUpperCase();
		return ("TAA".equals(value) || "TAG".equals(value) || "TGA".equals(value));
	}
}
