package org.vardb.util;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.codehaus.jackson.annotate.JsonProperty;

public class CCounts
{
	protected Map<String,List<Count>> types=new TreeMap<String,List<Count>>();

	public CCounts(){}
	
	public CCounts(List<CCountData> items)
	{
		for (CCountData item : items)
		{
			setCount(item);
		}
	}
	
	public void setCount(CCountData item)
	{
		setCount(item.getName(),item.getValue(),item.getCount().intValue());
	}
	
	public void setCount(String name, String value, int count)
	{
		//System.out.println("setting count: type="+type+", value="+value+", count="+count);
		findOrCreate(name,value).setCount(count);
	}
	
	@JsonProperty
	public Map<String,List<Count>> getTypes()
	{
		return types;
	}
	
	private Count findOrCreate(String name, String value)
	{
		List<Count> counts=findOrCreate(name);
		for (Count count : counts)
		{
			if (count.getValue().equals(value))
				return count;
		}
		Count count=new Count(name,value,0);
		counts.add(count);
		return count;
	}
	
	private List<Count> findOrCreate(String name)
	{
		if (!types.containsKey(name))
			types.put(name,new ArrayList<Count>());
		return types.get(name);
	}
	
	public class Count
	{
		protected String name;
		protected String value;		
		protected Integer count=null;
		
		public Count(String name, String value, int count)
		{
			this.name=name;
			this.value=value;
			this.count=count;
		}
		
		//@JsonProperty
		public String getName(){return this.name;}
		public void setName(final String name){this.name=name;}

		@JsonProperty
		public String getValue(){return this.value;}
		public void setValue(final String value){this.value=value;}

		@JsonProperty
		public Integer getCount(){return this.count;}
		public void setCount(final Integer count){this.count=count;}
	}
}
