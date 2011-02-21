package org.vardb.util.services;

import java.util.ArrayList;
import java.util.List;

import javax.servlet.http.HttpServletRequest;

import org.vardb.util.CFileHelper;
import org.vardb.util.CWebHelper;

public class FreemarkerHelper
{	
	protected String webapp;
	protected List<String> scriptpaths=new ArrayList<String>();
	protected List<String> scriptblocks=new ArrayList<String>();
	protected List<String> onreadylist=new ArrayList<String>();
	protected List<String> stylepaths=new ArrayList<String>();
	protected List<String> styleblocks=new ArrayList<String>();
	protected String timestamp=CFileHelper.getTimestamp();
	
	public FreemarkerHelper(HttpServletRequest request)
	{
		this.webapp=CWebHelper.getWebapp(request);
	}
	
	public void addRemoteScript(String url)
	{
		addRemoteScript(url,true);
	}
	
	public void addRemoteScript(String url, boolean addtimestamp)
	{
		if (addtimestamp)
			url=CWebHelper.appendParam(url, "param", timestamp);		
		this.scriptpaths.add(url);
	}
	
	public void addLocalScript(String filename)
	{
		addLocalScript(filename,true);
	}
	
	public void addLocalScript(String filename, boolean addtimestamp)
	{
		if (addtimestamp)
			filename=CWebHelper.appendParam(filename, "param", timestamp);
		this.scriptpaths.add(webapp+filename);
	}

	public void addScriptBlock(String block)
	{
		this.scriptblocks.add(block);
	}
	
	public void addOnReady(String onready)
	{
		onreadylist.add(onready);
	}
	
	public void onReady(String onready)
	{
		addOnReady(onready);
	}
	
	public String getScripts()
	{
		StringBuilder buffer=new StringBuilder();
		for (String path : scriptpaths)
		{
			buffer.append("<script src=\""+path+"\" type=\"text/javascript\"></script>\n");			
		}
		buffer.append("<script type=\"text/javascript\">\n");
		buffer.append("<!--\n");
		for (String block : scriptblocks)
		{
			buffer.append(block).append("\n");			
		}
		buffer.append("Ext.onReady(function()\n");
		buffer.append("{\n");
		for (String onready : onreadylist)
		{
			buffer.append(onready).append("\n");
		}
		buffer.append("});\n");
		buffer.append("// -->\n");
		buffer.append("</script>\n");
		return buffer.toString();
	}
	
	///////////////////////////////////////////////////
	
	public void addRemoteStylesheet(String url)
	{
		addRemoteStylesheet(url,true);
	}
	
	public void addRemoteStylesheet(String url, boolean addtimestamp)
	{
		if (addtimestamp)
			url=CWebHelper.appendParam(url, "param", timestamp);		
		this.stylepaths.add(url);
	}
	
	public void addLocalStylesheet(String filename)
	{
		addLocalStylesheet(filename,true);
	}
	
	public void addLocalStylesheet(String filename, boolean addtimestamp)
	{
		if (addtimestamp)
			filename=CWebHelper.appendParam(filename, "param", timestamp);
		this.stylepaths.add(webapp+filename);
	}

	public void addStyleBlock(String block)
	{
		this.styleblocks.add(block);
	}
	
	public String getCss()
	{
		StringBuilder buffer=new StringBuilder();
		for (String path : stylepaths)
		{
			buffer.append("<link rel=\"stylesheet\" type=\"text/css\" media=\"all\" href=\""+path+"\"/>\n");
		}
		buffer.append("<style type=\"text/css\">\n");
		for (String block : styleblocks)
		{
			buffer.append(block).append("\n");
		}
		buffer.append("</style>\n");
		return buffer.toString();
	}
}
