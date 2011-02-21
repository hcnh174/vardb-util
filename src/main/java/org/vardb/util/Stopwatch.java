package org.vardb.util;

// from http://www.goldb.org/stopwatchjava.html

public class Stopwatch
{
    private long startTime = 0;
    private long stopTime = 0;
    private boolean running = false;
    
    public Stopwatch(){}
    
    public Stopwatch(boolean start)
    {
    	if (start)
    		start();
    }
    
    public void start()
    {
        this.startTime = System.currentTimeMillis();
        this.running = true;
    }

    public long stop()
    {
        this.stopTime = System.currentTimeMillis();
        this.running = false;
        return this.getElapsedTime();
    }
    
    //elaspsed time in milliseconds
    public long getElapsedTime()
    {
        long elapsed;
        if (running)
        	elapsed = (System.currentTimeMillis() - startTime);
        else elapsed = (stopTime - startTime);
        return elapsed;
    }
    
    //elaspsed time in seconds
    public long getElapsedTimeSecs()
    {
        long elapsed;
        if (running)
        	elapsed = ((System.currentTimeMillis() - startTime) / 1000);
        else elapsed = ((stopTime - startTime) / 1000);
        return elapsed;
    }

    //sample usage
    public static void main(String[] args)
    {
        Stopwatch s = new Stopwatch();
        s.start();
        //code you want to time goes here
        s.stop();
       //System.out.println("elapsed time in milliseconds: " + s.getElapsedTime());
    }
}