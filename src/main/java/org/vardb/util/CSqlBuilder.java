package org.vardb.util;

import java.util.ArrayList;
import java.util.List;

public class CSqlBuilder
{
	protected List<String> withlist=new ArrayList<String>();
	protected List<String> selectlist=new ArrayList<String>();
	protected List<String> fromlist=new ArrayList<String>();
	protected List<String> aliaslist=new ArrayList<String>();
	protected List<String> wherelist=new ArrayList<String>();
	protected List<String> groupbylist=new ArrayList<String>();
	protected List<String> sortlist=new ArrayList<String>();
	protected Integer offset;
	protected Integer limit;

	public CSqlBuilder()
	{
		
	}
	
	public CSqlBuilder with(String value)
	{
		withlist.add(value);
		return this;
	}
	
	public CSqlBuilder select(String value)
	{
		selectlist.add(value);
		return this;
	}
	
	public CSqlBuilder from(String value)
	{
		fromlist.add(value);
		return this;
	}
	
	public CSqlBuilder where(String value)
	{
		wherelist.add(value);
		return this;
	}
	
	public CSqlBuilder groupby(String value)
	{
		groupbylist.add(value);
		return this;
	}
	
	public CSqlBuilder orderby(String value)
	{
		sortlist.add(value);
		return this;
	}
	
	public CSqlBuilder offset(int offset)
	{
		this.offset=offset;
		return this;
	}
	
	public CSqlBuilder limit(int limit)
	{
		this.limit=limit;
		return this;
	}
	
	@Override
	public String toString()
	{
		StringBuilder buffer=new StringBuilder();
		if (!this.withlist.isEmpty())
			buffer.append("with ").append(CStringHelper.join(withlist,",\n")).append("\n");
		buffer.append("select ").append(CStringHelper.join(selectlist,", ")).append("\n");
		buffer.append("from ").append(CStringHelper.join(fromlist,", ")).append("\n");
		if (!this.wherelist.isEmpty())
			buffer.append("where ").append(CStringHelper.join(wherelist," and ")).append("\n");
		if (!this.groupbylist.isEmpty())
			buffer.append("group by ").append(CStringHelper.join(groupbylist,", ")).append("\n");
		if (!this.sortlist.isEmpty())
			buffer.append("order by ").append(CStringHelper.join(sortlist,", ")).append("\n");
		if (this.offset!=null)
			buffer.append("offset "+this.offset+"\n");
		if (this.limit!=null) 
			buffer.append("limit "+this.limit+"\n");
		return buffer.toString();
	}
}