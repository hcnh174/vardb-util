package org.vardb.util;

import com.google.common.base.Throwables;

public final class CExceptionHelper
{	
	private CExceptionHelper(){}
	
	public static String getMessage(Exception e)
	{		
		StringBuilder buffer=new StringBuilder();
		Throwable t=(Throwable)e;
		while (t!=null)
		{
			String message=t.getMessage();
			if (message==null)
				message=t.toString();
			buffer.append(message);
			buffer.append("\n\n");
			t=t.getCause();
		}
		return buffer.toString();
	}
	
	public static String getFullStackTrace(Exception e)
	{
		StringBuilder buffer=new StringBuilder();
		Throwable t=(Throwable)e;
		while (t!=null)
		{
			getStackTrace(t,buffer);
			t=t.getCause();
		}
		return buffer.toString();
	}
	
	public static void getStackTrace(Throwable t, StringBuilder buffer)
	{
		StackTraceElement[] elements=t.getStackTrace();
		buffer.append(t.getClass().getName()+"\n");
		for (int index=0;index<elements.length;index++)
		{
			StackTraceElement element=elements[index];
			String classname=element.getClassName();
			String file=element.getFileName();
			String method=element.getMethodName();
			int line=element.getLineNumber();
	
			buffer.append("\t"+classname+"."+method+"("+file+":"+line+")\n");
		}
		buffer.append("\n\n");
	}
	
	public Throwable getRootCause(Throwable t)
	{
		return Throwables.getRootCause(t);
	}
}	