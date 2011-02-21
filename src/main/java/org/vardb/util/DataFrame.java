package org.vardb.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.vardb.util.CDataFrame.Column;

import com.google.common.collect.ForwardingTable;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;


/*
DataFrame<Integer,String,String> dataframe=new DataFrame<Integer,String,String>();//("patientid","name","birthdate");
dataframe.put(8888,"name","Bob Jones");
dataframe.put(8888,"birthdate","1/1/2001");
dataframe.put(9999,"name","Rita Smith");
dataframe.put(9999,"birthdate","9/9/1999");
System.out.println(dataframe.toString());
*/
public class DataFrame<R,C,V> extends ForwardingTable<R,C,V> {

	protected HashBasedTable<R,C,V> dataframe=HashBasedTable.create();
	protected List<String> colnames=new ArrayList<String>();
	
	
	public DataFrame(){}
	
	public DataFrame(String... colnames)
	{
		this.colnames.addAll(Arrays.asList(colnames));
	}
	
	public int getNumColumns()
	{
		return columnKeySet().size();
	}
	
	public int getNumRows()
	{
		return rowKeySet().size();
	}
	
	@Override
	public String toString()
	{
		StringBuilder buffer=new StringBuilder();
		boolean first=true;
		for (R rowname : rowKeySet())
		{
			Map<C,V> row=row(rowname);
			if (first)
			{
				buffer.append("id\t");
				buffer.append(CStringHelper.join(row.keySet(),"\t")).append("\n");
				first=false;
			}
			buffer.append(rowname).append("\t");
			buffer.append(CStringHelper.join(row.values(),"\t")).append("\n");
		}
		return buffer.toString();
	}
	
	/////////////////////////////////////////////
	
	@Override
	protected Table<R, C, V> delegate() {

		return dataframe;
	}

	
	
}
