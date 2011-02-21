package org.vardb.util;

import org.vardb.util.CDataType;

public class CAttributeList extends CValueList
{	
	public void addAttribute(String tagname, String name, String value)
	{
		addRow(tagname,name,value);
	}
	
	public void addAttributeType(String name, CDataType type, String description)
	{
		addRow(name,type.name(),description);
	}
}