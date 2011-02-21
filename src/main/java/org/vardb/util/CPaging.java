package org.vardb.util;

import org.vardb.util.CStringHelper;


public class CPaging implements IPaging
{
	protected boolean paged=false;
	protected Integer start;
	protected Integer pagesize;
	protected int total;
	//protected Filter filter=new Filter(); //hack!
	//protected CSorting sorting=new CSorting(); //hack!
	
	public boolean getPaged(){return this.paged;}
	public void setPaged(final boolean paged){this.paged=paged;}
	
	public Integer getStart(){return this.start;}
	public void setStart(Integer start){this.start=start;}

	public Integer getPagesize(){return this.pagesize;}
	public void setPagesize(Integer pagesize){this.pagesize=pagesize;}
	
	public int getTotal(){return this.total;}
	public void setTotal(int total){this.total=total;}
	
	//public ISorting getSorting(){return this.sorting;}
	//public IFilter getFilter(){return this.filter;}
	
	public CPaging()
	{
		this.start=0;
		this.pagesize=20;
		this.total=0;
	}
	
	public CPaging(Integer start, Integer pagesize)
	{
		this.start=(start==null) ? 0 : start;
		this.pagesize=(pagesize==null) ? 20 : pagesize;
	}
	
	//hack!
	public CPaging(Integer start, Integer pagesize, String filter, String sort, String dir)
	{
		this.start=(start==null) ? 0 : start;
		this.pagesize=(pagesize==null) ? 20 : pagesize;
		//this.filter.setFilter(filter);
		//this.sorting.setMultiSort(sort,dir);
	}
	
	@Override
	public String toString()
	{
		StringBuilder buffer=new StringBuilder();
		buffer.append("start="+this.start+"\n");
		buffer.append("pagesize="+this.pagesize+"\n");
		//buffer.append("sort="+this.sorting.getSql("")+"\n"); //hack!
		//buffer.append("filter="+this.filter+"\n"); //hack!
		return buffer.toString();
	}
	
	/*
	//hack!
	public class Filter implements IPaging.IFilter
	{
		protected String filter;
		
		public String getFilter(){return this.filter;}

		public Filter(){}
		
		public Filter(String filter)
		{
			this.filter=filter;
		}
		
		public void setFilter(String filter)
		{
			if (CStringHelper.isEmptyJson(filter))
				return;
			this.filter=filter;
		}
	}
	*/
}