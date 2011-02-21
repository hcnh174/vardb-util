package org.vardb.util;

import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.codehaus.jackson.annotate.JsonProperty;

import com.google.common.base.Charsets;
import com.google.common.base.Function;
import com.google.common.collect.BiMap;
import com.google.common.collect.Collections2;
import com.google.common.collect.HashBiMap;
import com.google.common.collect.Maps;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import com.google.common.io.LineReader;


//similar to a CTable but more geared to identifying cells using a double hash
// tables are row-based, dataframes are more column-based
public class CDataFrame
{
	protected Map<String,Column> columns=new LinkedHashMap<String,Column>();
	//protected Map<Object,Object> rownames=new LinkedHashMap<Object,Object>();
	protected Set<Object> rownames=new LinkedHashSet<Object>();
	protected boolean autoAddColumns=false;
	
	public CDataFrame(){}
	
	public CDataFrame(boolean autoAddColumns)
	{
		this.autoAddColumns=autoAddColumns;
	}
	
	public Column findOrCreateColumn(String colname)
	{
		Column column=getColumn(colname);
		if (column==null)
			column=addColumn(colname);
		return column;
	}
	
	public Column addColumn(String colname)
	{
		//System.out.println("adding column ["+colname+"]");
		if (this.columns.containsKey(colname))
			return this.columns.get(colname);
		Column column=new Column(this,colname);
		this.columns.put(column.getName(),column);
		return column;
	}

	public void addColumns(Collection<String> colnames)
	{
		for (String colname : colnames)
		{
			addColumn(colname);
		}
	}
	
	public Object getValue(String colname, Object rowname)
	{
		Column column=this.columns.get(colname);
		if (column==null)
			throw new CException("column is null: "+colname);
		return column.getValue(rowname);
	}
	
	public Object getValue(String colname, Object rowname, Object dflt)
	{
		Column column=this.columns.get(colname);
		if (column==null)
			return dflt;
		Object value=column.getValue(rowname);
		if (!CStringHelper.hasContent(value))
			return dflt;
		return value;
	}
	
	public void setValue(String colname, Object rowname, Object value)
	{
		//System.out.println("dataframe.setValue: colname="+colname+", rowname="+rowname+", value="+value);
		Column column=getColumn(colname);
		column.setValue(rowname,value);
	}
	
	public void setValue(String colname, Object rowname, char value)
	{
		Column column=getColumn(colname);
		column.setValue(rowname,value);
	}
	
	public boolean hasValue(String colname, Object rowname)
	{
		Object value=getValue(colname,rowname);
		return (value!=null);
	}
	
	public String getStringValue(String colname, Object rowname)
	{
		Object value=getValue(colname,rowname);
		if (value==null)
			return null;
		return value.toString();
	}
	
	public Integer getIntValue(String colname, Object rowname)
	{
		Object value=getValue(colname,rowname);
		if (value==null)
			return null;
		return Integer.valueOf(value.toString());
	}
	
	public Column getColumn(String colname)
	{
		if (!this.columns.containsKey(colname))
		{
			if (!autoAddColumns)
				throw new CException("can't find column named "+colname);
			else return addColumn(colname);
		}
		return this.columns.get(colname);
	}	
	
	public Collection<Column> getColumns()
	{
		return this.columns.values();
	}
	
	public int getNumRows()
	{
		return this.rownames.size();
	}
	
	public int getNumCols()
	{
		return this.columns.size();
	}
	
	public CTable getTable()
	{
		CTable table=new CTable();
		//table.getHeader().add("rowname");
		for (String colname : getColNames())
		{
			table.getHeader().add(colname);
		}
		for (Object rowname : getRowNames())
		{
			CTable.Row row=table.addRow();
			//row.add(rowname);
			for (String colname : getColNames())
			{
				Column column=getColumn(colname);
				Object value=column.getValue(rowname);
				row.add(value);
			}
		}
		return table;
	}

	public List<String> getColNames()
	{
		List<String> colnames=new ArrayList<String>();
		colnames.addAll(this.columns.keySet());
		return colnames;
	}
	
	public boolean hasRow(Object rowname)
	{
		return this.rownames.contains(rowname);
		//return this.rownames.containsKey(rowname);
	}
	
	public boolean hasColumn(String colname)
	{
		return this.columns.containsKey(colname);
	}
	
	public void registerRowName(Object rowname)
	{
		addRow(rowname);
	}
	
	public void addRow(Object rowname)
	{
		if (!hasRow(rowname))
			this.rownames.add(rowname);
			//this.rownames.put(rowname,rowname);
	}
	
	public Collection<Object> getRowNames()
	{
		return rownames;
		//return rownames.keySet();
	}
	
	public Collection<String> getStringRowNames()
	{
		//return Collections2.transform(rownames.keySet(), new Function<Object, String>()
		return Collections2.transform(rownames, new Function<Object, String>()
		{
			@Override
			public String apply(Object obj)
			{
				return obj.toString();
			}
		});
	}
	
	public int size()
	{
		Column column=this.columns.values().iterator().next();
		return column.size();
	}
	
	
	public boolean isEmpty()
	{
		return this.columns.isEmpty();
	}
	
	@Override
	public String toString()
	{
		return getTable().toString();
	}
	
	public void appendColumns(CDataFrame other)
	{
		for (String colname : other.getColNames())
		{
			if (hasColumn(colname))
			{
				System.err.println("dataframe already has column "+colname+". skipping.");
				continue;
			}
			addColumn(colname);
			CDataFrame.Column column=other.getColumn(colname);
			for (Object rowname : column.getRowNames())
			{
				Object value=column.getValue(rowname);
				if (hasRow(rowname))
					setValue(colname,rowname,value);
			}
		}
	}
	
	public void appendForeignColumns(String fkey_colname, CDataFrame other)
	{
		for (String other_colname : other.getColNames())
		{
			String colname=fkey_colname+"_"+other_colname;
			if (hasColumn(colname))
			{
				System.err.println("dataframe already has column "+colname+". skipping.");
				continue;
			}
			addColumn(colname);
			CDataFrame.Column column=other.getColumn(other_colname);
			for (Object other_rowname : column.getRowNames())
			{
				Object value=column.getValue(other_rowname);
				for (Object rowname : this.getRownamesByColValue(fkey_colname, other_rowname))
				{
					//System.out.println("trying to setValue(colname="+colname+", rowname="+rowname+", value="+value);
					setValue(colname,rowname,value);
				}
			}
		}
	}
	
	public Map<String,CDataType> getColTypes()
	{
		Map<String,CDataType> coltypes=new LinkedHashMap<String,CDataType>();
		for (CDataFrame.Column column : getColumns())
		{
			coltypes.put(column.getName(),column.guessDataType());
		}
		return coltypes;
	}
	
	public Collection<Object> getUniqueValues(String colname)
	{
		CDataFrame.Column column=getColumn(colname);
		return column.getUniqueValues();
	}
	
	public Collection<Object> getRownamesByColValue(String colname, Object value)
	{
		CDataFrame.Column column=getColumn(colname);
		return column.getRownamesByColValue(value); 
	}
	
	public CIdList getUniqueIds(String colname)
	{
		CDataFrame.Column column=getColumn(colname);
		return column.getUniqueIds();
	}

	public CAttributeList getAttributeList()
	{
		CAttributeList attlist=new CAttributeList();
		for (Object tagname : getRowNames())
		{
			for (String attname : getColNames())
			{
				Object value=getValue(attname,tagname);
				attlist.addAttribute(tagname.toString(), attname, value.toString());
			}
		}
		return attlist;
	}
	
	@JsonProperty
	public Integer getTotalCount()
	{
		return rownames.size();
	}
	
	@JsonProperty
	public List<Row> getRows()
	{
		List<Row> rows=new ArrayList<Row>();
		for (Object rowname : getRowNames())
		{
			Row row=new Row();
			row.put("rowname",rowname);
			for (String colname : getColNames())
			{
				Object value=getValue(colname,rowname);
				row.put(colname,value);
			}
			rows.add(row);
		}
		return rows;
	}
	
	@SuppressWarnings("serial")
	public static class Row extends LinkedHashMap<String,Object>
	{
		
	}
	
	public static class Column
	{
		protected CDataFrame dataframe;
		protected String colname;
		protected Map<Object,Object> values=new LinkedHashMap<Object,Object>();
		
		public Column(CDataFrame dataframe, String colname)
		{
			this.dataframe=dataframe;
			this.colname=colname;
		}
		
		public String getName(){return this.colname;}
		
		public Set<Object> getKeys()
		{
			return this.values.keySet();
		}
		
		public Object getValue(Object rowname)
		{
			return this.values.get(rowname);
		}
		
		public void setValue(Object rowname, Object value)
		{
			dataframe.registerRowName(rowname);
			this.values.put(rowname,value);
		}
		
		public Collection<Object> getRowNames()
		{
			return this.values.keySet();
		}
		
		public Collection<Object> getValues()
		{
			return this.values.values();
		}
		
		public int size()
		{
			return this.values.size();
		}
		
		public CDataType guessDataType()
		{
			return CDataType.guessDataType(values.values());
		}
		
		public Collection<Object> getUniqueValues()
		{
			Set<Object> uniquevalues=new LinkedHashSet<Object>();
			for (Object value : values.values())
			{
				if (!uniquevalues.contains(value))
					uniquevalues.add(value);
			}
			return uniquevalues;
		}

		public CIdList getUniqueIds()
		{
			Set<Integer> uniquevalues=new LinkedHashSet<Integer>();
			for (Object val : values.values())
			{
				Integer value=Integer.valueOf(val.toString());
				if (!uniquevalues.contains(value))
					uniquevalues.add(value);
			}
			return new CIdList(uniquevalues);
		}
		
		public Collection<Object> getRownamesByColValue(Object val)
		{
			Set<Object> rownames=new LinkedHashSet<Object>();
			for (Object rowname : values.keySet())
			{
				Object value=values.get(rowname);
				System.out.println("if "+value+".equals("+val+"): "+value.equals(val));
				if (value.equals(val))
					rownames.add(rowname);
			}
			return rownames;
		}
	}

	public interface IntervalListener
	{
		public boolean onInterval(Parser parser, int rownum);
	}
	
	public abstract static class Parser
	{
		protected final int ROW_STATUS_INTERVAL=100000;
		protected int interval=10;
		protected Charset encoding=Charsets.UTF_8;
		protected List<IntervalListener> listeners=new ArrayList<IntervalListener>();
		protected BiMap<Integer,String> headerfields=HashBiMap.create();
		protected CDataFrame dataframe;
		protected Map<Integer,CDataFrame.Column> columns;//=new LinkedHashMap<Integer,CDataFrame.Column>();
		
		public Parser(){}
		
		public Parser(Charset encoding)
		{
			this.encoding=encoding;
		}
		
		public int getInterval(){return this.interval;}
		public void setInterval(final int interval){this.interval=interval;}
		
		public void addIntervalListener(IntervalListener listener)
		{
			listeners.add(listener);
		}
		
		public CDataFrame getDataFrame(){return dataframe;}
		
		protected boolean readHeader(List<String> fields)
		{
			fields=preProcessHeader(fields);
			for (int index=0;index<fields.size();index++)
			{
				String colname=fields.get(index);
				headerfields.put(index, colname);
			}
			postProcessHeader(fields);
			resetDataFrame();
			return true;
		}
		
		public void resetDataFrame()
		{
			dataframe=new CDataFrame();
			columns=Maps.newLinkedHashMap();
			for (Integer index : headerfields.keySet())
			{
				String colname=headerfields.get(index);
				Column column=this.dataframe.addColumn(colname);
				this.columns.put(index,column);
			}
		}
		
		/*
		protected boolean readHeader(List<String> fields)
		{
			fields=preProcessHeader(fields);
			for (int index=0;index<fields.size();index++)
			{
				String colname=fields.get(index);
				Column column=this.dataframe.addColumn(colname);
				this.columns.put(index,column);
			}
			postProcessHeader(fields);
			return true;
		}
		*/
		protected boolean readLine(List<String> values, int rownum)
		{
			values=preProcessLine(values,rownum);
			if (values.size()!=columns.size())
				throw new CException("numbers of fields and headings don't match: fields="+values.size()+", columns="+columns.size()+" in row "+rownum);
			String rowname=values.get(0).trim(); // assume the first column is the row name
			//this.dataframe.rownames.put(rowname,rowname);
			this.dataframe.rownames.add(rowname);
			for (int index=0;index<values.size();index++)
			{
				Column column=this.columns.get(index);
				String value=values.get(index).trim();
				//System.out.println("setting col="+column.getName()+", row="+rownum+", value="+value);
				column.values.put(rowname,value);
			}
			postProcessLine(values,rownum);
			return true;
		}
		
		protected boolean notifyListeners(int rownum)
		{
			if (rownum%interval!=0)
				return true;
			boolean proceed=true;
			for (IntervalListener listener : listeners)
			{
				if (!listener.onInterval(this, rownum))
					proceed=false;
			}
			return proceed;
		}
		
		protected List<String> preProcessHeader(List<String> fields){return fields;}
		protected void postProcessHeader(List<String> fields){}
		
		protected List<String> preProcessLine(List<String> values, int rownum)
		{
			return values;
		}
	
		protected List<String> splitLine(String line)
		{
			return CStringHelper.splitAsList(line,"\t");//line.split("\t");
		}
		
		protected void postProcessLine(List<String> values, int rownum)
		{
			if (rownum%ROW_STATUS_INTERVAL==0)
				System.out.println("reading line "+rownum);
		}
	}
	
	public static class TabFileParser extends Parser
	{	
		public TabFileParser(){}
		
		public TabFileParser(Charset encoding)
		{
			super(encoding);
		}
		
		public void parseFile(String filename)
		{
			try		
			{
				Files.readLines(new File(filename), this.encoding, new LineProcessor<String>()
				{
					private int rownum=0;
					
					public boolean processLine(String line)
					{
						//System.out.println("line: "+line);
						if (rownum==0)
							readHeader(splitLine(line));
						else readLine(splitLine(line),rownum);
						if (!notifyListeners(rownum))
							return false;
						rownum++;
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
		
		public void parse(String str)
		{
			try		
			{
				LineReader reader=new LineReader(new StringReader(str));
				String line;
				int rownum=0;
				while ((line=reader.readLine())!=null)
				{
					System.out.println("reading line: "+line);
					if (rownum==0)
						readHeader(splitLine(line));
					else readLine(splitLine(line),rownum);
					rownum++;
				}
			}
			catch (IOException e)
			{
				throw new CException(e);
			}
		}
	}
	
	public static CDataFrame parseTabFile(String filename)
	{
		return parseTabFile(filename,Charsets.UTF_8);
	}
	
	public static CDataFrame parseTabFile(String filename, Charset encoding)
	{
		TabFileParser parser=new TabFileParser(encoding);
		parser.parseFile(filename);
		return parser.getDataFrame();
	}
	
	public static CDataFrame parse(String str)
	{
		return parse(str,Charsets.UTF_8);
	}
	
	public static CDataFrame parse(String str, Charset encoding)
	{
		TabFileParser parser=new TabFileParser(encoding);
		parser.parse(str);
		return parser.getDataFrame();
	}
}