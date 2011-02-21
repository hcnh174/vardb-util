package org.vardb.util;

import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.Locale;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public final class CDateHelper
{
	//public static final String DATE_PATTERN="MM'/'dd'/'yyyy";
	//public static final String DATE_PATTERN="yyyy'/'dd'/'MM";
	//public static final String DATE_PATTERN="MM'/'dd'/'yyyy";
	//public static final String DATETIME_PATTERN=DATE_PATTERN+"' '"+TIME_PATTERN;
	
	public static final String MMDDYYYY_PATTERN="MM'/'dd'/'yyyy";
	public static final String YYYYMMDD_PATTERN="yyyy'/'MM'/'dd";
	public static final String POSTGRES_YYYYMMDD_PATTERN="yyyy'-'MM'-'dd";	//										  
	public static final String EXCEL_PATTERN="EEE MMM dd hh:mm:ss zzz yyyy";//Thu Oct 07 00:00:00 JST 1982
	public static final String DATE_PATTERN=MMDDYYYY_PATTERN;
	public static final String TIME_PATTERN="h':'mm' 'a";
	public static final String DATETIME_PATTERN=DATE_PATTERN+"' 'H':'mm':'ss";
	private static final Locale LOCALE=Locale.US;
	
	private CDateHelper(){}
	
	public static String format(String pattern)
	{
		return format(new Date(),pattern);
	}
	
	// Formatting
	public static String format(Date date, String pattern)
	{
		if (date==null)
			return "";
		SimpleDateFormat formatter=new SimpleDateFormat(pattern,LOCALE);
		return formatter.format(date);
	}
	
	public static String format(Date date)
	{
		//String pattern=DATE_PATTERN+"' '"+TIME_PATTERN;
		String pattern=DATETIME_PATTERN;
		return format(date,pattern);
	}
	
	public static Date parse(String value)
	{
		return parse(value,DATETIME_PATTERN);
	}
	
	public static Date parse(String value, String pattern)
	{
		return parse(value,pattern,true);
	}
	
	// try to handle some common data cleaning problems
	public static Date parse(String origvalue, String pattern, boolean throwException)
	{
		try
		{
			//String value=cleanDate(origvalue);
			String value=origvalue;
			if (value==null)
				return null;
			//if (value.indexOf(" ")!=-1)
			//	pattern=DATETIME_PATTERN; // EXCEL_PATTERN //hack!
			SimpleDateFormat formatter=new SimpleDateFormat(pattern,LOCALE);
			return formatter.parse(value.trim());
		}
		catch (ParseException e)
		{
			if (throwException)
				throw new CException("problem parsing date ["+origvalue+"] with pattern ["+pattern+"]",e);
			if (CStringHelper.hasContent(origvalue))
				System.err.println("problem parsing date ["+origvalue+"] with pattern ["+pattern+"]");
			return null;
		}
	}
	
	private static final String YEAR_REGEX="([0-9]{4})?";
	private static final String YEAR_MONTH_REGEX="([0-9]{4})/([0-9]{1,2})?";
	private static final String FLEXIBLE_DATE_REGEX="([0-9]{4})[^0-9]([0-9]+)[^0-9]([0-9]+)";
	
	public static String cleanDate(String origvalue)
	{
		String value=origvalue;
		if ("??".equals(value) || "????".equals(value))
			return null;
		
		// remove extraneous characters
		if (value.indexOf('?')!=-1)
			value=value.replace("?","");
		if (value.indexOf("//")!=-1)
			value=value.replace("//","/"); // 2005//6/13
				
		// try to determine pattern
		if (value.matches(YEAR_REGEX))
			return value.substring(0,4)+"/1/1";
		else if (value.matches(YEAR_MONTH_REGEX)) //1994/10?
		{	
			Pattern pat=Pattern.compile(YEAR_MONTH_REGEX);
			Matcher matcher=pat.matcher(value);
			matcher.find();
			String year=String.valueOf(Integer.parseInt(matcher.group(1)));
			String month=String.valueOf(Integer.parseInt(matcher.group(2)));
			return year+"/"+month+"/1";
		}
		else if (value.matches(FLEXIBLE_DATE_REGEX))
		{
			Pattern pat=Pattern.compile(FLEXIBLE_DATE_REGEX);
			Matcher matcher=pat.matcher(value);
			matcher.find();
			String year=String.valueOf(Integer.parseInt(matcher.group(1)));
			String month=String.valueOf(Integer.parseInt(matcher.group(2)));
			String date=String.valueOf(Integer.parseInt(matcher.group(3)));
			return year+"/"+month+"/"+date;
		}
		return value;
	}
	
	public static boolean isDate(String value, String pattern)
	{
		try
		{
			SimpleDateFormat formatter=new SimpleDateFormat(pattern,LOCALE);
			formatter.parse(value);
			return true;
		}
		catch (ParseException e)
		{
			return false;
		}
	}
	
	public static Date setDate(String strdate)
	{
		return parse(strdate,DATE_PATTERN);
	}

	public static Date setDate(String strdate, String strtime)
	{
		return parse(strdate+"T"+strtime,DATE_PATTERN+"'T'"+TIME_PATTERN);
	}
	
	public static Date setDate(int year, int month, int date)
	{
		Calendar cal=new GregorianCalendar();
		cal.set(Calendar.YEAR,year);
		cal.set(Calendar.MONTH,month);
		cal.set(Calendar.DATE,date);
		return cal.getTime();
	}
	
	public static Date setTime(Date date, String strtime)
	{
		Date time=parse(strtime,TIME_PATTERN);
		return setTime(date,time);
	}
	
	public static Date setTime(String strtime)
	{
		Date date=new Date();
		return setTime(date,strtime);
	}
	
	public static Date setTime(Date datepart, Date timepart)
	{
		Calendar date=new GregorianCalendar();
		Calendar time=new GregorianCalendar();
		date.setTime(datepart);
		time.setTime(timepart);
		date.set(Calendar.HOUR_OF_DAY,time.get(Calendar.HOUR_OF_DAY));
		date.set(Calendar.MINUTE,time.get(Calendar.MINUTE));
		date.set(Calendar.SECOND,time.get(Calendar.SECOND));
		date.set(Calendar.MILLISECOND,0);
		return date.getTime();
	}
	
	public static Date addDays(Date date, int days)
	{
		Calendar calendar=new GregorianCalendar();
		calendar.setTime(date);
		calendar.add(Calendar.DATE,days);
		return calendar.getTime();
	}
	
	public static Date addMinutes(Date date, int minutes)
	{
		Calendar calendar=new GregorianCalendar();
		calendar.setTime(date);
		calendar.add(Calendar.MINUTE,minutes);
		return calendar.getTime();
	}
	
	public static Date addHours(Date date, int hours)
	{
		Calendar calendar=new GregorianCalendar();
		calendar.setTime(date);
		calendar.add(Calendar.HOUR,hours);
		return calendar.getTime();
	}
	
	public static Date addWeeks(Date date, int weeks)
	{ 
		return addDays(date,weeks*7);
	}
	
	public static Date addMonths(Date date, int months)
	{
		Calendar calendar=new GregorianCalendar();
		calendar.setTime(date);
		calendar.add(Calendar.MONTH,months);
		return calendar.getTime();
	}
	
	public static Date setToMidnight(Date date)
	{
		return setToHour(date,0);
	}
	
	public static Date setToNearestHour(Date date)
	{
		Calendar cal=new GregorianCalendar();
		cal.setTime(date);
		return setToHour(date,cal.get(Calendar.HOUR_OF_DAY));
	}
	
	public static Date setToHour(Date date, int hour)
	{
		Calendar cal=new GregorianCalendar();
		cal.setTime(date);		
		Calendar calendar=new GregorianCalendar();
		calendar.set(cal.get(Calendar.YEAR),cal.get(Calendar.MONTH),cal.get(Calendar.DATE));
		calendar.set(Calendar.HOUR_OF_DAY,hour);
		calendar.set(Calendar.MINUTE,0);
		calendar.set(Calendar.SECOND,0);
		calendar.set(Calendar.MILLISECOND,0);
		return calendar.getTime();
	}
	
	/*
	public static Date setToMidnight(Date date)
	{
		Calendar calendar=new GregorianCalendar();
		calendar.setTime(date);
		calendar.set(Calendar.HOUR_OF_DAY,0);
		calendar.set(Calendar.MINUTE,0);
		calendar.set(Calendar.SECOND,0);
		return calendar.getTime();
	}
	*/
	
	public static Date stripDate(Date date)
	{
		Calendar calendar=new GregorianCalendar();
		calendar.setTime(date);
		calendar.set(Calendar.YEAR,1900);
		calendar.set(Calendar.MONTH,0);
		calendar.set(Calendar.DATE,1);
		return calendar.getTime();
	}
	
	public static Date getWeekStartDate(Date date, int startday)
	{
		Calendar calendar=new GregorianCalendar();
		calendar.setTime(date);
		int curday=calendar.get(Calendar.DAY_OF_WEEK);
		int days=curday-startday;
		days=-1*days;
		return addDays(date,days);
	}
	
	public static Date getFirstDayInMonth(Date date)
	{
		Calendar calendar=new GregorianCalendar();
		calendar.setTime(date);
		calendar.set(Calendar.DATE,1);
		return calendar.getTime();
	}
	
	public static int getYear(Date date)
	{
		Calendar calendar=new GregorianCalendar();
		calendar.setTime(date);
		return calendar.get(Calendar.YEAR);
	}
	
	public static int getMonth(Date date)
	{
		Calendar calendar=new GregorianCalendar();
		calendar.setTime(date);
		return calendar.get(Calendar.MONTH);
	}
	
	public static int getDate(Date date)
	{
		Calendar calendar=new GregorianCalendar();
		calendar.setTime(date);
		return calendar.get(Calendar.DAY_OF_MONTH);
	}
	
	public static int getDay(Date date)
	{
		Calendar calendar=new GregorianCalendar();
		calendar.setTime(date);
		return calendar.get(Calendar.DAY_OF_WEEK);	
	}
	
	// returns duration in minutes
	public static int getDuration(Date date1, Date date2)
	{
		long diff=date2.getTime()-date1.getTime();
		return (int)(diff/(60*1000));
	}
	
	public static Integer getDurationInYears(Date date1, Date date2)
	{
		if (date1==null || date2==null)
		{
			System.err.println("can't determine duration because of null date. date1="+date1+", date2="+date2);
			return null;
		}
		Calendar cal1=new GregorianCalendar();
		Calendar cal2=new GregorianCalendar();
		cal1.setTime(date1);
		cal2.setTime(date2);
		int year1=cal1.get(Calendar.YEAR);
		int year2=cal2.get(Calendar.YEAR);
		return Math.abs(year2-year1);		
	}
	
	public static int getDaysInMonth(Date date)
	{
		Calendar calendar=new GregorianCalendar();
		calendar.setTime(date);
		calendar.set(Calendar.DATE,1);
		calendar.add(Calendar.MONTH,1);
		calendar.add(Calendar.DATE,-1);
		return calendar.get(Calendar.DATE);
	}
	
	// if startdate is null, sets start to midnight on current date and enddate to one day from then
	public static class DateRange
	{
		protected Date startdate;
		protected Date enddate;
		
		public Date getStartdate(){return this.startdate;}
		public Date getEnddate(){return this.enddate;}
		
		public DateRange(Date startdate, Date enddate)
		{
			this.startdate=startdate;
			this.enddate=enddate;
		}
		
		public DateRange(String strstartdate, String strenddate)
		{
			this(strstartdate,strenddate,CDateHelper.DATE_PATTERN);
		}
		 
		public DateRange(String strstartdate, String strenddate, String pattern)
		{
			if (strstartdate==null)
			{
				this.startdate=CDateHelper.setToMidnight(new Date());
				this.enddate=CDateHelper.addDays(this.startdate,1);
			}
			else
			{
				this.startdate=CDateHelper.parse(strstartdate,pattern);
				this.enddate=CDateHelper.parse(strenddate,pattern);
			}
		}
	}
}