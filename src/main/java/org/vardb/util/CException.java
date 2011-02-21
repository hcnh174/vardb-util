package org.vardb.util;

@SuppressWarnings("serial")
public class CException extends RuntimeException
{
	public CException()
	{
		super();
	}
	
	public CException(Throwable t)
	{
		super(t);
	}
	
	public CException(String message)
	{
		super(message);
	}
	
	public CException(String message, Throwable t)
	{
		super(message,t);
	}
}