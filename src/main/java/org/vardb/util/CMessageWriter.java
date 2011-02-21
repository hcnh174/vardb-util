package org.vardb.util;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.List;

import javax.servlet.http.HttpServletResponse;

import org.vardb.util.CException;
import org.vardb.util.CStringHelper;

public class CMessageWriter
{
	protected PrintWriter writer;
	protected int errors=0;
	protected boolean isSystemOut=false;
	
	public CMessageWriter(HttpServletResponse response)
	{
		try
		{
			this.writer=response.getWriter();
			isSystemOut=isSystemOut();
		}
		catch(Exception e)
		{
			throw new CException(e);
		}
	}
	
	public CMessageWriter(PrintWriter writer)
	{
		this.writer=writer;
		isSystemOut=isSystemOut();
	}
	
	public CMessageWriter()
	{
		this.writer=new PrintWriter(System.out);
		isSystemOut=true;
	}
	
	public void setQuiet()
	{
		writer=new PrintWriter(new StringWriter());
	}
	
	public void write(String str)
	{
		//logger.info(str);
		//if (!isSystemOut)
		//	System.out.print(str);
		this.writer.print(str);
	}
	
	public void br()
	{
		if (!isSystemOut)
			this.writer.print("<br>\n");
		else this.writer.println("");
	}
	
	public void message(String message)
	{
		if (!isSystemOut)
			System.out.println(message);
		write(message);
		br();
		flush();
	}
	
	public void messages(List<String> messages)
	{
		for (String message : messages)
		{
			message(message);
		}
	}
	
	public void warn(String message)
	{
		this.errors++;
		message("<span style=\"color:red;\">WARN: "+message+"</span>");
	}
	
	public void error(String message)
	{
		this.errors++;
		message("<span style=\"color:red;\">ERROR: "+message+"</span>");
	}
	
	public void error(Exception e)
	{
		this.errors++;
		error(e.toString());
		CStringHelper.println(e.toString());
		e.printStackTrace();
	}
	
	public void error(String message, Exception e)
	{
		error(message);
		error(e);
	}
	
	public PrintWriter getWriter(){return this.writer;}
	
	public void flush()
	{
		this.writer.flush();
	}
	
	public boolean hasErrors()
	{
		return this.errors>0;
	}
	
	
	//////////////////////////////////////////
	
	private boolean isSystemOut()
	{
		System.out.println("writer: "+writer.getClass().getCanonicalName());
		System.out.println("system.out: "+System.out.getClass().getCanonicalName());
		return writer.getClass().getCanonicalName().equals(System.out.getClass().getCanonicalName());
	}
}