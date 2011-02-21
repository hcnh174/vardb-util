package org.vardb.util;

import java.io.ByteArrayInputStream;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Map;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.vardb.util.CXmlValidationException.ValidationError;
import org.xml.sax.ErrorHandler;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.SAXParseException;

public final class CXmlHelper
{
	private CXmlHelper(){}
	
	public static void start(StringBuilder buffer, String name)
	{
		buffer.append("<"+name);
	}
	
	public static void end(StringBuilder buffer, String name)
	{
		buffer.append("</"+name+">");
	}
	
	public static void element(StringBuilder buffer, String name, String value)
	{
		if (CStringHelper.hasContent(value))
			buffer.append("\t<"+name+">"+value+"</"+name+">\n");
	}
	
	public static void cdata(StringBuilder buffer, String name, String value)
	{
		if (CStringHelper.hasContent(value))
			buffer.append("\t<"+name+"><![CDATA["+value+"]]></"+name+">\n");
	}
	
	public static void attribute(StringBuilder buffer, String name, Object value)
	{
		if (CStringHelper.hasContent(value))
			buffer.append(" "+name+"=\""+value.toString()+"\"");
	}
	
	public static void attribute(StringBuilder buffer, String name, Object value, Object dlft)
	{
		if (!CStringHelper.hasContent(value))
			value=dlft;
		attribute(buffer,name,value);
	}
	
	public static void attributes(StringBuilder buffer, Map<String,String> attributes)
	{
		for (String name : attributes.keySet())
		{
			String value=attributes.get(name);
			attribute(buffer,name,value);
		}
	}

	private static final String JAXP_SCHEMA_LANGUAGE = "http://java.sun.com/xml/jaxp/properties/schemaLanguage";
	private static final String JAXP_SCHEMA_SOURCE = "http://java.sun.com/xml/jaxp/properties/schemaSource";	
	private static final String W3C_XML_SCHEMA = "http://www.w3.org/2001/XMLSchema";
	
	public static void validate(String xml, String xsd)	//String xsdfile)
		throws CXmlValidationException
	{
		CErrorHandler handler=new CErrorHandler();
		try
		{
			DocumentBuilderFactory factory=DocumentBuilderFactory.newInstance();
			factory.setNamespaceAware(true);
			factory.setValidating(true);
		
			factory.setAttribute(JAXP_SCHEMA_LANGUAGE, W3C_XML_SCHEMA);
			factory.setAttribute(JAXP_SCHEMA_SOURCE, new ByteArrayInputStream(xsd.getBytes()));
			//factory.setAttribute(JAXP_SCHEMA_SOURCE, new File(xsdfile));

			DocumentBuilder parser=factory.newDocumentBuilder();
			parser.setErrorHandler(handler);
			parser.parse(new InputSource(new StringReader(xml)));
		}
		catch (Exception e)
		{
			throw new CException(e);
		}
		if (handler.hasErrors())
			handler.throwException();
	}
	
	public static class CErrorHandler implements ErrorHandler
	{
		protected List<ValidationError> m_errors=new ArrayList<ValidationError>();
		
		public void warning(SAXParseException e)
			throws SAXException
		{
			add(e);
		}
	    
	    public void error(SAXParseException e)
	    	throws SAXException
	    {
	    	add(e);
	    }
	    
	    public void fatalError(SAXParseException e)
	    	throws SAXException
	    {
	    	add(e);
	    }
	    
	    ////////////////////////////////////////////
	   
	    private void add(SAXParseException e)
	    {
	    	m_errors.add(new ValidationError(e));
	    }
	    
	    public boolean hasErrors()
	    {
	    	return !m_errors.isEmpty();
	    }
	    
	    public void throwException() throws CXmlValidationException
	    {
	    	throw new CXmlValidationException(m_errors);
	    }
	}
	
	///////////////////////////////////////////////////////////////////////
	
	public static String removePI(String xml, String element)
	{
		int index=xml.indexOf("<"+element);
		if (index!=-1)
			xml=xml.substring(index);
		return xml.trim();
	}
	
	public static String removeRootElement(String xml, String element)
	{
		int index=xml.indexOf("<"+element+">");
		if (index!=-1)
		{
			xml=xml.substring(index); // remove xml declaration, doctype, etc.
			xml=xml.replace("<"+element+">","");
			xml=xml.replace("</"+element+">","");
		}
		return xml.trim();
	}
	
	public static String addRootElement(String xml, String element)
	{
		return "<"+element+">"+xml+"</"+element+">";
	}
	
	public static String removeDoctype(String xml)
	{
		int start=xml.indexOf("<!DOCTYPE");
		if (start==-1)
			return xml;
		int end=xml.indexOf('>',start)+1;
		xml=xml.substring(0,start)+xml.substring(end+1);
		if (xml.indexOf("<!DOCTYPE")!=-1)
			throw new CException("DOCTYPE still in xml: "+xml);
		return xml;
	}
	
	////////////////////////////////////////////////////////////////
	
	public static String mergeXmlFiles(String folder, String root)
	{
		List<String> filenames=CFileHelper.listFilesRecursively(folder,".xml");
		return mergeXmlFiles(filenames,root);
	}
	
	public static String mergeXmlFiles(String folder, String root, Date date)
	{
		List<String> filenames=CFileHelper.listFilesRecursively(folder,".xml",date);
		return mergeXmlFiles(filenames,root);
	}
	
	public static String mergeXmlFiles(List<String> filenames, String root)
	{
		StringBuilder buffer=new StringBuilder();
		buffer.append("<"+root+">\n");
		for (String filename : filenames)
		{
			String xml=CFileHelper.readFile(filename);
			xml=removeRootElement(xml,root);
			buffer.append("\n\n");
			buffer.append(xml);
		}
		buffer.append("</"+root+">\n");
		CFileHelper.writeFile("c:/temp/setup.xml",buffer.toString());
		return buffer.toString();
	}
	
	// validate
	public static void main(String[] argv)
	{
		String dir=argv[0];
		String xsdfile=argv[1];
		String xsd=CFileHelper.readFile(xsdfile);
		for (String filename : CFileHelper.listFilesRecursively(dir,".xml"))
		{
			try
			{
				//System.out.println("filename="+filename);
				String xml=CFileHelper.readFile(filename);
				CXmlHelper.validate(xml,xsd);
			}
			catch(CXmlValidationException e)
			{
				System.out.println(filename);
				System.out.println(e.getMessage());
			}
			catch(Exception e)
			{
				System.out.println(filename);
				System.out.println(e.getMessage());
				e.printStackTrace();
			}
		}
	}
}