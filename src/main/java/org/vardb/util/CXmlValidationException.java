package org.vardb.util;

import java.util.ArrayList;
import java.util.List;

import org.xml.sax.SAXParseException;

@SuppressWarnings("serial")
public class CXmlValidationException extends RuntimeException
{
	protected List<ValidationError> errors=new ArrayList<ValidationError>();
	
	public CXmlValidationException(List<ValidationError> list)
	{
		super();
		add(list);
	}
	
	public List<ValidationError> getErrors(){return this.errors;}
	
	@Override
	public String getMessage()
	{
		StringBuilder buffer=new StringBuilder();
		buffer.append("Xml validation errors\n");
		for (ValidationError error : this.errors)
		{
			buffer.append(error.toString()).append("\n");
		}
		return buffer.toString();
	}
	
	private void add(List<ValidationError> list)
    {
		for (ValidationError e : list)
		{
			add(e);
		}
    }
	
	private void add(ValidationError e)
    {
    	this.errors.add(e);
    }
	
	public static class ValidationError
	{
		protected Integer lineNumber;
		protected String message;
		
		public Integer getLineNumber(){return this.lineNumber;}
		public String getMessage(){return this.message;}

		public ValidationError(SAXParseException e)
		{
			this.lineNumber=e.getLineNumber();
			this.message=e.getMessage();
		}
		
		@Override
		public String toString()
		{
			return "Error on line "+this.lineNumber+": "+this.message;
		}
	}
}