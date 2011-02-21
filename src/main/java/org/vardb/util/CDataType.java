package org.vardb.util;

import java.util.Collection;
import java.util.Date;

public enum CDataType
{
	STRING(false,"string"),
	INTEGER(true,"int"),
	FLOAT(true,"float"),
	BOOLEAN(false,"boolean"),
	DATE(false,"date");
	
	private boolean numeric;
	private String json;
	
	CDataType(boolean numeric, String json)
	{
		this.numeric=numeric;
		this.json=json;
	}
	
	public boolean isNumeric(){return this.numeric;}
	public String getJson(){return this.json;}
	
	public static CDataType guessDataType(Collection<Object> values)
	{
		boolean is_integer=true;
		boolean is_float=true;
		boolean is_boolean=true;			
		for (Object obj : values)
		{
			String value=obj.toString().toLowerCase().trim();
			if (CStringHelper.isEmpty(value))
				continue;
			if (!CMathHelper.isFloat(value))
				is_float=false;
			if (!CMathHelper.isInteger(value))
				is_integer=false;
			if (!"true".equals(value) && !"false".equals(value))
				is_boolean=false;
		}
		if (is_boolean)
			return CDataType.BOOLEAN;
		if (is_integer)
			return CDataType.INTEGER;
		if (is_float)
			return CDataType.FLOAT;			
		return CDataType.STRING;
	}
	
	public static CDataType guessDataType(Object obj)
	{
		CDataType type=guessDataTypeByClass(obj);
		if (obj==null || type!=CDataType.STRING)
			return type;
		boolean is_integer=true;
		boolean is_float=true;
		boolean is_boolean=true;
		String value=obj.toString().toLowerCase().trim();
		if (CStringHelper.isEmpty(value))
			return CDataType.STRING;
		if (!CMathHelper.isFloat(value))
			is_float=false;
		if (!CMathHelper.isInteger(value))
			is_integer=false;
		if (!"true".equals(value) && !"false".equals(value))
			is_boolean=false;
		// determine data type
		if (is_boolean)
			return CDataType.BOOLEAN;
		if (is_integer)
			return CDataType.INTEGER;
		if (is_float)
			return CDataType.FLOAT;			
		return CDataType.STRING;
	}
	
	private static CDataType guessDataTypeByClass(Object obj)
	{
		if (obj instanceof String)
			return CDataType.STRING;
		else if (obj instanceof Integer)
			return CDataType.INTEGER;
		else if (obj instanceof Float || obj instanceof Double)
			return CDataType.FLOAT;
		else if (obj instanceof Boolean)
			return CDataType.BOOLEAN;
		else if (obj instanceof Date)
			return CDataType.DATE;
		else if (obj instanceof Byte)
			return CDataType.STRING;
		if (obj!=null)
			throw new CException("no guessDataType handler for class "+obj.getClass().getName()+": ["+obj.toString()+"]");
			//System.err.println("no guessDataType handler for class "+obj.getClass().getName());
		return CDataType.STRING;
	}
}
