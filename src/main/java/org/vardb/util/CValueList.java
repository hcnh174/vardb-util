package org.vardb.util;

import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import org.vardb.util.CStringHelper;

import com.google.common.base.Function;
import com.google.common.collect.Collections2;
import com.google.common.collect.Lists;

public class CValueList
{
	protected List<Row> rows=Lists.newArrayList();
	
	public void addRow(String... args)
	{
		addRow(new Row(args));
	}
	
	public void addRow(Row row)
	{
		rows.add(row);
	}
	
	public boolean isEmpty()
	{
		return rows.isEmpty();
	}
	
	public int size()
	{
		return rows.size();
	}
	
	public String getSql(String alias)
	{
		StringBuilder buffer=new StringBuilder();
		buffer.append("(\n");
		buffer.append("values\n");
		buffer.append(CStringHelper.join(rows,",\n"));
		buffer.append("\n) as "+alias+"\n");
		return buffer.toString();
	}
	
	// creates a new attributelist based on unique values from the first column only
	public CAttributeList getUnique()
	{
		CAttributeList list=new CAttributeList();
		Set<String> values=new LinkedHashSet<String>();
		for (Row row : rows)
		{
			values.add(row.columns.get(0));
		}
		for (String value : values)
		{
			list.addRow(value);
		}
		return list;
	}
	
	public class Row
	{
		protected List<String> columns;
		
		public Row(String... args)
		{
			columns=Arrays.asList(args);
		}
		
		@Override
		public String toString()
		{
			Collection<String> values=Collections2.transform(columns,new Function<String,String>()
			{
				@Override
				public String apply(String value)
				{
					return CStringHelper.singleQuote(CStringHelper.escapeSql(value));
				}
				
			});
			StringBuilder buffer=new StringBuilder();
			buffer.append("(");
			buffer.append(CStringHelper.join(values,","));
			buffer.append(")");
			return buffer.toString();
		}
	}
}