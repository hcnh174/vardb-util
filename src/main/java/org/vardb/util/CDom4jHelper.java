package org.vardb.util;

import java.io.StringReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.dom4j.Attribute;
import org.dom4j.Document;
import org.dom4j.Element;
import org.dom4j.Node;
import org.dom4j.io.SAXReader;

import com.google.common.collect.Sets;

public final class CDom4jHelper
{	
	private static final String ELEMENT_IS_NULL="element is null";
	
	private CDom4jHelper(){}
	
	public static Document parse(String xml)
	{
		try
		{
			SAXReader reader = new SAXReader();
			return reader.read(new StringReader(xml));
		}
		catch(Exception  e)
		{
			throw new CException(e);
		}
	}
	
	public static String getAttribute(Element element, String name)
	{
		if (element==null)
			throw new CException(ELEMENT_IS_NULL);
		Attribute attribute=element.attribute(name);
		if (attribute==null)
			throw new CException("attribute "+name+" is null");
		String value=attribute.getValue();
		//System.out.println("attribute="+value);
		return fixEntities(value);
	}
	
	// optional default value
	public static String getAttribute(Element element, String name, Object dflt)
	{
		if (element==null)
			throw new CException(ELEMENT_IS_NULL);
		Attribute attribute=element.attribute(name);
		if (attribute==null)
		{
			if (dflt==null)
				return null;
			else return dflt.toString();
		}
		String value=attribute.getValue();
		//System.out.println("attribute="+value);
		return fixEntities(value);
	}
	
	public static Integer getIntAttribute(Element element, String name)
	{
		String value=getAttribute(element,name);
		return Integer.valueOf(value);
	}
	
	public static Integer getIntAttribute(Element element, String name, Integer dflt)
	{
		String value=getAttribute(element,name,dflt);
		return Integer.valueOf(value);
	}
	
	public static Boolean getBoolAttribute(Element element, String name)
	{
		String value=getAttribute(element,name);
		return Boolean.valueOf(value);
	}
	
	public static Boolean getBoolAttribute(Element element, String name, Boolean dflt)
	{
		String value=getAttribute(element,name,dflt);
		return Boolean.valueOf(value);
	}

	public static String getValue(Element element, String path)
	{
		return getValue(element,path,null);
	}
	
	public static String getValue(Element element, String path, String dflt)
	{
		if (element==null)
			throw new CException(ELEMENT_IS_NULL);
		String value=element.valueOf(path);
		if (value==null)
		{
			if (dflt==null)
				throw new CException("value is null for path"+path);
			else return dflt;
		}
		//System.out.println("path value="+value);
		return fixEntities(value);
	}
	
	public static String getText(Node node)
	{
		//System.out.println("node.getText()="+node.getText());
		if (node==null)
			throw new CException("node is null");
		String text=node.getText();
		return fixEntities(text);
	}
	
	public static String getTrimmedText(Node node)
	{
		return getText(node).trim();
	}
	
	public static String getChildrenAsXml(Element element)
	{
		if (element==null)
			throw new CException(ELEMENT_IS_NULL);
		StringBuilder buffer=new StringBuilder();
		Iterator<?> iter=element.nodeIterator();
		while(iter.hasNext())
		{
			Node child = (Node)iter.next();
			//System.out.println("child name="+child.getName());
			//System.out.println("string value=["+child.getStringValue().trim()+"]");
			//System.out.println("xml=["+child.asXML()+"]");
			//buffer.append(child.asXML().trim());
			buffer.append(child.asXML());
		}
		String value=buffer.toString().trim();
		//String value=buffer.toString();
		return fixEntities(value);
	}
	
	public static String getChildrenAsXml(Element parent, String path)
	{
		Element element=(Element)parent.selectSingleNode(path);
		return getChildrenAsXml(element);
	}
	
	public static String getChildrenAsXml(Element parent, String path, String dflt)
	{
		Element element=(Element)parent.selectSingleNode(path);
		if (element==null)
			return dflt;
		return getChildrenAsXml(element);
	}
		
	public static String fixEntities(String str)
	{
		//System.out.println("trying to fix entities="+str);
		//return CStringHelper.replace(str,AMPERSAND,"&");
		return str;
	}

	/*
	public static List<Element> selectNodes(Element root, String name)
	{
		List<Element> nodes=new ArrayList<Element>();
		for (Iterator<?> iter=root.nodeIterator();iter.hasNext();)
		{
			Object node=iter.next();
			if (!(node instanceof Element))
				continue;
			Element element=(Element)node;
			if (!element.getName().equals(name))
				continue;
			nodes.add(element);
		}
		System.out.println("found "+nodes.size()+" nodes of type"+name);
		return nodes;
	}
	*/
	
	public static List<Element> selectNodes(Element root, String... args)
	{
		Set<String> names=Sets.newHashSet(args);
		List<Element> nodes=new ArrayList<Element>();
		for (Iterator<?> iter=root.nodeIterator();iter.hasNext();)
		{
			Object node=iter.next();
			if (!(node instanceof Element))
				continue;
			Element element=(Element)node;
			if (!names.contains(element.getName()))
				continue;
			nodes.add(element);
		}
		System.out.println("found "+nodes.size()+" nodes of type"+names);
		return nodes;
	}
}