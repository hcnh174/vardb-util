package org.vardb.util;


public final class ThreadHelper
{	
	private ThreadHelper(){}
	
	public static void sleep(long millis)
	{
		System.out.println("going to sleep: "+millis+" milliseconds");
		Stopwatch stopwatch=new Stopwatch();
		stopwatch.start();
		try
		{
			Thread.sleep(millis);
		}
		catch (InterruptedException e)
		{
			CStringHelper.println("sleep interrupted: "+e);
		}
		System.out.println("waking up after "+stopwatch.stop()+" milliseconds");
	}
}