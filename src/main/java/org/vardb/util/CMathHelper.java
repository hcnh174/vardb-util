package org.vardb.util;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public final class CMathHelper
{
	private CMathHelper(){}
	
	public static int randomInteger(int size)
	{
		//int rnd=(int)(double)(Math.random()*(double)size);
		//return rnd;
		Random random=new Random();
		return random.nextInt(size);
	}
	
	public static double calculateDiversity(List<Integer> counts)
	{
		int total=0;
		for (Integer count : counts)
		{
			total+=count;
		}
		
		if (total==0)
			return 0;
		
		double sum=0.0f;
		for (Integer count : counts)
		{
			double p=((double)count/(double)total);
			sum+=p*p;
		}
		return 1/sum;
	}
	
	public static boolean isInteger(String str)
	{
		if (str==null)
			return false;
		str=str.trim();
		if (CStringHelper.isEmpty(str))
			return false;
		final char[] numbers = str.toCharArray();
		for (int x = 0; x < numbers.length; x++)
		{      
			final char c = numbers[x];
			if ((c >= '0') && (c <= '9'))
				continue;
			return false; // invalid
		}
		return true; // valid
	}
	
	/*
	public static boolean isFloat(String str)
	{
		if (str==null)
			return false;
		str=str.trim();
		if (CStringHelper.isEmpty(str))
			return false;
		final char[] numbers = str.toCharArray();
		int numdots=0;
		for (int x = 0; x < numbers.length; x++)
		{      
			final char c = numbers[x];
			if ((c >= '0') && (c <= '9'))
				continue;
			if (c=='.')
			{
				numdots++;
				continue;
			}
			return false; // invalid
		}
		if (numdots>1) // if there is more than one decimal point
			return false;
		return true; // valid
	}
	*/
	
	public static Integer parseInt(String str, int dflt)
	{
		Integer num=parseInt(str);
		if (num==null)
			return dflt;
		return num;
	}
	
	public static Integer parseInt(String str)
	{
		if (!CStringHelper.hasContent(str))
			return null;
		try
		{
			if (str.matches("[0-9]+\\.0+")) // remove empty decimal places
				str=str.substring(0,str.indexOf('.'));
			return Integer.valueOf(str);
		}
		catch(NumberFormatException e)
		{
			return null;
		}
	}
	
	public static boolean isFloat(String str)
	{
		return parseFloat(str)!=null;
	}
	
	public static Float parseFloat(String str)
	{
		//if (!CStringHelper.hasContent(str))
		//	return null;
		try
		{
			return Float.valueOf(str);
		}
		catch(NumberFormatException e)
		{
			return null;
		}
	}
	
	public static Double parseDouble(String str)
	{
		if (!CStringHelper.hasContent(str))
			return null;
		try
		{
			return Double.valueOf(str);
		}
		catch(NumberFormatException e)
		{
			return null;
		}
	}
	
	// tries to parse a boolean field from text
	// considers + equivalent to true and - equivalent to false
	public static Boolean parseBoolean(String str)
	{
		if (!CStringHelper.hasContent(str))
			return null;
		if ("+".equals(str))
			return true;
		else if ("-".equals(str))
			return false;
		try
		{
			return Boolean.valueOf(str);
		}
		catch(NumberFormatException e)
		{
			return null;
		}
	}
	
	public static List<Integer> getIntersection(List<List<Integer>> lists)
	{
		List<Integer> ids=new ArrayList<Integer>();
		for (List<Integer> list : lists)
		{
			for (Integer id : list)
			{
				if (ids.contains(id))
					continue;
				if (intersects(lists, id))
					ids.add(id);
			}
		}
		return ids;
	}
	
	private static boolean intersects(List<List<Integer>> lists, int id)
	{
		if (lists.isEmpty())
			return false;
		boolean intersects=true;
		for (List<Integer> list : lists)
		{
			if (!list.contains(id))
				intersects=false;
		}
		return intersects;
	}
	
	public static double log2(double num)
	{
		if (num==0)
			throw new CException("can't take log of 0");
		return (Math.log(num)/Math.log(2));
	}
	
	public static int getNumbatches(int total, int batchsize)
	{
		int numbatches=(int)Math.floor((double)total/(double)batchsize);
		if (((double)total)%((double)batchsize)!=0)
			numbatches++;
		System.out.println("numbatches="+numbatches);
		return numbatches;
	}
	
	public static List<List<String>> getBatches(List<String> ids, int batchsize)
	{
		int iterations=(int)Math.floor((double)ids.size()/(double)batchsize);
		List<List<String>> batches=new ArrayList<List<String>>();
		int start=0;
		int end=0;
		for (int batchnumber=0;batchnumber<iterations;batchnumber++)
		{
			start=batchnumber*batchsize;
			end=start+batchsize;
			//System.out.println("batch - from "+start+" to "+end);
			List<String> batch=ids.subList(start,end);
			batches.add(batch);
		}
		if (((double)ids.size())%((double)batchsize)!=0)
			batches.add(ids.subList(start,ids.size()));
		return batches;
	}
}
