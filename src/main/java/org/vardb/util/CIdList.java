package org.vardb.util
;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import javax.persistence.MappedSuperclass;
import javax.persistence.Transient;

@MappedSuperclass
public class CIdList
{
	private static final int MIN_RANGE_SIZE=10;
	private static final int MAX_LIST_SIZE=500;
	private static final String RANGE_DELIMITER="-";
	private static final String LIST_DELIMITER=",";
	
	private String idlist="";
	private List<Range> ranges=new ArrayList<Range>();

	public String getIdlist(){return this.idlist;}
	public void setIdlist(final String idlist){this.idlist=idlist;}
	
	public CIdList() {}

	public CIdList(int id)
	{
		addId(id);
	}
	
	public CIdList(Collection<Integer> ids)
	{
		addIds(ids);
	}
	
	public CIdList(CIdList ids)
	{
		addIds(ids);
	}
	
	public CIdList(String str)
	{
		parse(str);
	}

	public void onSave()
	{
		clean();
		if (!this.ranges.isEmpty())
			this.idlist=toString();
		else this.idlist="";
		//System.out.println("CSerializableIdList.onSave called: ranges="+this.ranges.size()+" ids="+this.idlist);
	}
	
	public void onLoad()
	{
		parse(this.idlist);
		//System.out.println("CSerializableIdList.onLoad called: ranges="+this.ranges.size()+" ids="+this.idlist);
	}
	
	private void clean()
	{
		Collections.sort(this.ranges,new RangeComparator());
	}
	
	private void parse(String idlist)
	{
		this.ranges.clear();
		if (!CStringHelper.hasContent(idlist))
			return;
		for (String str : CStringHelper.split(idlist,LIST_DELIMITER))
		{
			addRange(new Range(str));
		}
		Collections.sort(this.ranges,new RangeComparator());
	}
	
	public void clear()
	{
		this.ranges.clear();
		this.idlist="";
	}
	
	@Transient
	public boolean isEmpty()
	{
		return this.ranges.isEmpty();
	}
	
	@Transient
	public int getSize()
	{
		return size();
	}
	
	public int size()
	{
		int count=0;
		for (Range range : this.ranges)
		{
			count+=range.size();
		}
		return count;
	}
	
	/*
	@Transient
	public CResourceIDs getResourceIds(CResourceType type)
	{
		return new CResourceIDs(type,new CIdList(toString()));
	}
	*/
	
	public void setIds(List<Integer> ids)
	{
		clear();
		addIds(ids);
	}
	
	public void setIds(CIdList ids)
	{
		clear();
		addIds(ids);
	}

	public boolean contains(int id)
	{
		for (Range range : this.ranges)
		{
			if (range.overlaps(id))
				return true;
		}
		return false;
	}
	
	public int addIds(CIdList ids)
	{
		return addIds(ids.getIds());
	}
	
	public int addIds(Collection<Integer> rawids)
	{
		//System.out.println("CIdList.addIds: "+rawids.size()+" ids");
		if (rawids.isEmpty())
			return 0;
		List<Integer> ids=asList(rawids);
		List<Range> ranges=new ArrayList<Range>();
		List<Integer> other=new ArrayList<Integer>();		
		findRanges(ids,ranges,other);
		int numadded=0;
		for (Range range : ranges)
		{
			numadded+=addRange(range);
		}		
		for (Integer id : other)
		{
			numadded+=addId(id);
		}
		mergeRanges();
		return numadded;
	}
	
	private void findRanges(List<Integer> ids, List<Range> ranges, List<Integer> other)
	{
		Collections.sort(ids);	
		// hack - add a -1 to start of list to simplify code
		ids.add(0,-1);
		ids.add(Integer.MAX_VALUE);
		
		List<List<Integer>> lists=new ArrayList<List<Integer>>();
		List<Integer> list=new ArrayList<Integer>();		
		for (int index=1;index<ids.size()-1;index++)
		{
			int lastid=ids.get(index-1);
			int id=ids.get(index);
			int nextid=ids.get(index+1);
			if (lastid==id-1 && nextid==id+1) // part of a series
			{
				list.add(id);
			}
			else if (lastid!=id-1 && nextid==id+1) // start a new series
			{
				list=new ArrayList<Integer>();
				list.add(id);
			}
			else if (lastid==id-1 && nextid!=id+1) // end of a series
			{
				list.add(id);
				lists.add(list);
			}
			else if (lastid!=id-1 && nextid!=id+1) // end of a series
			{
				other.add(id);
			}
		}
		for (List<Integer> sublist : lists)
		{
			ranges.add(new Range(sublist));
		}
	}

	private int addRange(Range range)
	{
		for (Range r : this.ranges)
		{
			if (r.overlaps(range))
				return r.merge(range);
		}
		this.ranges.add(range);
		return range.size();
	}
	
	private void mergeRanges()
	{
		Collections.sort(this.ranges,new RangeComparator());
		List<Range> ranges=new ArrayList<Range>();
		Range range=this.ranges.get(0);
		ranges.add(range);
		for (int index=1;index<this.ranges.size();index++)
		{
			Range next=this.ranges.get(index);
			if (range.overlaps(next))
			{
				range.merge(next);
				continue;
			}
			if (range.adjacent(next))
			{
				range.extend(next);
				continue;
			}
			ranges.add(next);
			range=next;
		}
		this.ranges=ranges;
	}
	
	public int addId(int id)
	{
		for (Range r : this.ranges)
		{
			if (r.overlaps(id))
				return 0;
			if (r.adjacent(id))
			{
				r.extend(id);
				return 1;
			}
		}
		addRange(new Range(id,id));
		return 1;
	}
	
	public void deleteId(int id)
	{
		deleteIds(Collections.singletonList(id));
	}
	
	public int deleteIds(CIdList ids)
	{
		return deleteIds(ids.getIds());
	}
	
	public int deleteIds(Collection<Integer> deletelist)
	{
		//System.out.println("CIdList.deleteIds: "+deletelist.size()+" ids");
		List<Integer> ids=new ArrayList<Integer>();
		for (Integer id : getIds())
		{
			boolean keep=true;
			for (Integer deleteid : deletelist)
			{
				int id1=id;
				int id2=deleteid;
				if (id1==id2)
				{
					keep=false;
					break;
				}
			}
			if (keep)
				ids.add(id);
		}
		Collections.sort(ids);
		int numremoved=getIds().size()-ids.size();
		setIds(ids);
		return numremoved;
	}
	
	@Transient
	public List<Integer> getIds()
	{
		List<Integer> list=new ArrayList<Integer>();
		for (Range range : this.ranges)
		{
			range.getIds(list);
		}
		return list;
	}
	
	@Transient
	private List<Integer> getDiscontinuous()
	{
		List<Integer> list=new ArrayList<Integer>();
		for (Range range : this.ranges)
		{
			if (range.size()<MIN_RANGE_SIZE)
				range.getIds(list);
		}
		return list;
	}
	
	public String toSql()
	{
		return toSql("sequence.id");
	}
	
	public String toSql(String field)
	{
		clean();
		if (size()==0)
			return "(1=2)";
		List<String> subqueries=new ArrayList<String>();
		for (Range range : this.ranges)
		{
			if (range.size()>=MIN_RANGE_SIZE)
				addSubquery(subqueries,range.toSql(field));
		}
		
		List<Integer> other=getDiscontinuous();		
		createInListSubquery(field,other,subqueries);// put the remaining ids in id list, splitting if necessary
		String sql=CStringHelper.join(subqueries," or ").trim();
		return CStringHelper.parenthesize(sql);
	}
	
	private void createInListSubquery(String field, Collection<Integer> ids, List<String> subqueries)
	{
		// split the ids into smaller groups
		Collection<Collection<Integer>> lists=CStringHelper.split(ids,MAX_LIST_SIZE);
		for (Collection<Integer> list : lists)
		{
			addSubquery(subqueries,createInList(field,list));
		}
	}
	
	private String createInList(String field, Collection<Integer> ids)
	{
		if (ids.isEmpty())
			return "";
		if (ids.size()==1)
			return field+"="+ids.iterator().next();
		else return field+" in ("+CStringHelper.join(ids,",")+")";
	}
	
	///////////////////////////////////////////////////
	
	// convert it to a unique list
	private List<Integer> asList(Collection<Integer> col)
	{
		List<Integer> list=new ArrayList<Integer>();
		for (Integer id : col)
		{
			if (!list.contains(id))
				list.add(id);
		}
		return list;
	}
	
	private static void addSubquery(List<String> subqueries, String subquery)
	{
		if (CStringHelper.hasContent(subquery))
			subqueries.add(subquery);
	}
	
	public boolean intersects(CIdList other)
	{
		for (int id : other.getIds())
		{
			if (contains(id))
				return true;
		}
		return false;
	}
	
	public boolean identical(CIdList other)
	{
		return toString().equals(other.toString());
	}
	
	@Override
	public String toString()
	{
		StringBuilder buffer=new StringBuilder();
		toString(buffer);
		return buffer.toString();
	}
	
	public void toString(StringBuilder buffer)
	{
		clean();
		if (this.ranges.isEmpty())
		{
			buffer.append("");
			return;
		}
		for (int index=0;index<this.ranges.size();index++)
		{
			if (index!=0)
				buffer.append(LIST_DELIMITER);
			Range range=this.ranges.get(index);
			range.toString(buffer);
		}
	}
	
	public String toLocation()
	{
		StringBuilder buffer=new StringBuilder();
		for (int index=0;index<this.ranges.size();index++)
		{
			if (index!=0)
				buffer.append(LIST_DELIMITER);
			Range range=this.ranges.get(index);
			range.toLocation(buffer);
		}
		return buffer.toString();
	}
	
	public static void initialize(Collection<? extends CIdList> idlists)
	{
		for (CIdList idlist : idlists)
		{
			idlist.onLoad();
		}
	}
	
	public static class Range
	{
		private int lower;
		private int upper;
		
		public int getLower(){return this.lower;}
		public void setLower(final int lower){this.lower=lower;}

		public int getUpper(){return this.upper;}
		public void setUpper(final int upper){this.upper=upper;}
		
		public Range(int lower, int upper)
		{
			this.lower=lower;
			this.upper=upper;
			check();
		}
		
		public Range(int id)
		{
			this(id,id);
		}		
		
		public Range(List<Integer> ids)
		{
			this(ids.get(0),ids.get(ids.size()-1));
		}
		
		public Range(String str)
		{
			int index=str.indexOf(RANGE_DELIMITER);
			if (index==-1)
			{
				int id=Integer.parseInt(str);
				this.lower=id;
				this.upper=id;
				return;
			}
			this.lower=Integer.parseInt(str.substring(0,index));
			this.upper=Integer.parseInt(str.substring(index+1));
			check();
		}
		
		public void check()
		{
			if (this.lower>this.upper)
			{
				//System.out.println("swapping range: lower="+this.lower+", upper="+this.upper);
				int temp=this.upper;
				this.upper=this.lower;
				this.lower=temp;
			}
		}
		
		public boolean overlaps(Range other)
		{
			return (overlaps(other.getLower()) || overlaps(other.getUpper()));
		}
		
		public boolean overlaps(int id)
		{
			return (id>=this.lower && id <=this.upper);
		}
		
		public boolean adjacent(Range other)
		{
			return ((this.upper+1==other.getLower()) || (other.getUpper()+1==this.lower));
		}
		
		public boolean adjacent(int id)
		{
			return (id==this.lower-1 || id==this.upper+1);
			//return (id==this.upper+1);
		}
		
		public void extend(Range other)
		{
			if (this.upper+1==other.getLower())
				this.upper=other.getUpper();
			else if (other.getUpper()+1==this.lower)
				this.lower=other.getLower();
		}
		
		public void extend(int id)
		{
			if (id==this.lower-1)
				this.lower=id;
			else if (id==this.upper+1)
				this.upper=id;
		}
		
		public int merge(Range other)
		{
			this.lower=(this.lower<other.getLower()) ? this.lower : other.getLower();
			this.upper=(this.upper>other.getUpper()) ? this.upper : other.getUpper();
			return intersection(other).size();
		}
		
		public Range intersection(Range other)
		{
			// find whichever lower value is greater
			int lower=(this.lower>other.getLower()) ? this.lower : other.getLower();
			int upper=(this.upper<other.getUpper()) ? this.upper : other.getUpper();
			return new Range(lower,upper);
		}
		
		
		public void split(int id, List<Range> ranges)
		{
			ranges.add(new Range(this.lower,id-1));
			ranges.add(new Range(id+1,this.upper));
		}
		
		public void split(Range range, List<Range> ranges)
		{
			Range range1=new Range(this.lower,range.getLower()-1);
			Range range2=new Range(range.getUpper()+1,this.upper);
			//System.out.println("split, new range1: "+range1.toString()+", range2: "+range2.toString());
			ranges.add(range1);
			ranges.add(range2);
		}
		
		public int size()
		{
			return this.upper-this.lower+1;
		}
		
		public boolean isSingleton()
		{
			return (this.lower==this.upper);
		}
		
		public void getIds(List<Integer> ids)
		{
			for (int id=this.lower;id<=this.upper;id++)
			{
				ids.add(id);
			}
		}
		
		public String toSql(String field)
		{
			return field+" between "+this.lower+" and "+this.upper;
		}
		
		public String toString()
		{
			StringBuilder buffer=new StringBuilder();
			toString(buffer);
			return buffer.toString();
		}
		
		public void toString(StringBuilder buffer)
		{
			if (isSingleton())
				buffer.append(this.lower);
			else buffer.append(this.lower).append(RANGE_DELIMITER).append(this.upper);
		}
		
		// the only differences are that the delimiter should be ".." and singletons should be like 4..4
		public void toLocation(StringBuilder buffer)
		{
			if (isSingleton())
				buffer.append(this.lower+".."+this.lower);
			else buffer.append(this.lower).append("..").append(this.upper);
		}
	}

	@SuppressWarnings("serial")
	public static class RangeComparator implements Comparator<Range>, Serializable
	{
		public int compare(Range r1, Range r2)
		{
			//return r2.getLower()-r1.getLower();
			return r1.getLower()-r2.getLower();
		}
	}
	
	
}