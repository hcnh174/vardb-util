package org.vardb.util;

import java.beans.PropertyDescriptor;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import org.springframework.beans.BeanUtils;
import org.springframework.beans.BeanWrapper;
import org.springframework.beans.ConfigurablePropertyAccessor;
import org.springframework.beans.PropertyAccessorFactory;
import org.springframework.beans.propertyeditors.CustomDateEditor;

public class CBeanHelper
{
	protected String logfile="c:/temp/beanutilerrors.txt";
	protected boolean logerrors=true;
	
	public void copyProperties(Object target, Object src)
	{
		if (target==null)
			throw new CException("CBeanHelper.copy(): target object is null");
		if (src==null)
			throw new CException("CBeanHelper.copy(): source object is null");
		try
		{
			BeanUtils.copyProperties(src, target);
		}
		catch(Exception e)
		{
			CFileHelper.appendFile(logfile,e.getMessage());
		}
	}
	
	public Object getProperty(Object target, String property)
	{
		try
		{
			BeanWrapper wrapper=PropertyAccessorFactory.forBeanPropertyAccess(target);
			return wrapper.getPropertyValue(property);
		}
		catch(Exception e)
		{
			CFileHelper.appendFile(logfile,e.getMessage());
			return null;
		}
	}
	
	public boolean setPropertyFromString(Object target, String property, String value)
	{
		try
		{
			value=CStringHelper.trim(value);
			BeanWrapper wrapper=PropertyAccessorFactory.forBeanPropertyAccess(target);
			wrapper.setPropertyValue(property,value);
			return true;
		}
		catch(Exception e)
		{
			CFileHelper.appendFile(logfile,e.getMessage());
			return false;
		}
	}
	
	public boolean setPropertyFromString(Object target, String property, String value, String datepattern)
	{
		try
		{
			value=CStringHelper.trim(value);
			BeanWrapper wrapper=PropertyAccessorFactory.forBeanPropertyAccess(target);
			wrapper.registerCustomEditor(Date.class,new CustomDateEditor(new SimpleDateFormat(datepattern),true));
			wrapper.setPropertyValue(property,value);
			return true;
		}
		catch(Exception e)
		{
			CFileHelper.appendFile(logfile,e.getMessage());
			return false;
		}
	}
	
	public static List<String> getProperties(Object obj)
	{
		List<String> properties=new ArrayList<String>();
		try
		{
			BeanWrapper wrapper=PropertyAccessorFactory.forBeanPropertyAccess(obj);
			for (PropertyDescriptor property : wrapper.getPropertyDescriptors())
			{
				properties.add(property.getName());
			}
			return properties;
		}
		catch(Exception e)
		{
			throw new CException(e);
			//CFileHelper.appendFile(logfile,e.getMessage());
		}
	}
	
	public static Object forName(String cls)
	{
		try
		{
			return Class.forName(cls);
		}
		catch (ClassNotFoundException e)
		{
			throw new CException(e);
		}
	}
	
	@SuppressWarnings("unchecked")
	public static Object newInstance(Class cls)
	{
		try
		{
			return cls.newInstance();
		}
		catch(IllegalAccessException e)
		{
			throw new CException(e);
		}
		catch(InstantiationException e)
		{
			throw new CException(e);
		}
	}
}