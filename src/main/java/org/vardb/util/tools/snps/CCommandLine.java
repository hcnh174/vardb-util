package org.vardb.util.tools.snps;

import java.io.File;
import java.util.List;

import org.apache.commons.exec.CommandLine;
import org.apache.commons.exec.DefaultExecutor;
import org.apache.commons.exec.ExecuteWatchdog;
import org.vardb.util.CPlatformType;
import org.vardb.util.CStringHelper;

public class CCommandLine
{
	protected CommandLine commandLine=null;
	protected File workingDir=null;
	protected Integer watchdog=null;
	protected Integer exitValue=null;
	protected boolean throwExceptionIfOutputStream=false;
	protected boolean throwExceptionIfErrorStream=false;
	protected boolean cygwin=false;
	protected CPlatformType.Platform platform=CPlatformType.find().getPlatform();

	public CCommandLine(String command)
	{
		this.commandLine = CommandLine.parse(command);
	}
	
	public CCommandLine(List<String> commands)
	{
		addArgs(commands);
	}

	public CommandLine getCommandLine()
	{
		if (this.cygwin && this.platform==CPlatformType.Platform.WINDOWS)
			return getCygwinCommandLine(this.commandLine);
		return this.commandLine;
	}
	
	public static CommandLine getCygwinCommandLine(CommandLine commandLine)
	{
		CCommandLine commands=new CCommandLine("bash");
		commands.addArg("--login");
		commands.addArg("-c");
		commands.addArg(commandLine.toString(),true);
		return commands.getCommandLine();
	}
	
	public void addArg(Object value)
	{
		if (value==null)
			return;
		this.commandLine.addArgument(value.toString());
	}
	
	public void addArg(Object value, boolean handleQuoting)
	{
		if (value==null)
			return;
		this.commandLine.addArgument(value.toString(),handleQuoting);
	}
	
	public void addArgIf(boolean addarg, Object value)
	{
		if (!addarg)
			return;
		if (value==null)
			return;
		this.commandLine.addArgument(value.toString());
	}
	
	public void addArgs(List<String> args)
	{
		for (String arg : args)
		{
			addArg(arg);
		}
	}
	
	public void addArg(String name, Object value)
	{
		if (value==null)
			return;
		addArg(name);
		addArg(value);
	}
	
	public void addArgIf(boolean addarg, String name, Object value)
	{
		if (!addarg)
			return;
		if (value==null)
			return;
		addArg(name);
		addArg(value);
	}
	
	public void setWorkingDir(String dir)
	{
		setWorkingDir(new File(dir));
	}
	
	public File getWorkingDir(){return this.workingDir;}
	public void setWorkingDir(final File workingDir){this.workingDir=workingDir;}
	
	public Integer getWatchdog(){return this.watchdog;}
	public void setWatchdog(int watchdog){this.watchdog=watchdog;}
	
	public Integer getExitValue(){return this.exitValue;}
	public void setExitValue(int exitValue){this.exitValue=exitValue;}
	
	public boolean getThrowExceptionIfOutputStream(){return this.throwExceptionIfOutputStream;}
	public void setThrowExceptionIfOutputStream(final boolean throwExceptionIfOutputStream){this.throwExceptionIfOutputStream=throwExceptionIfOutputStream;}
	
	public boolean getThrowExceptionIfErrorStream(){return this.throwExceptionIfErrorStream;}
	public void setThrowExceptionIfErrorStream(final boolean throwExceptionIfErrorStream){this.throwExceptionIfErrorStream=throwExceptionIfErrorStream;}
	
	public boolean getCygwin(){return this.cygwin;}
	public void setCygwin(final boolean cygwin){this.cygwin=cygwin;}
	
	public CPlatformType.Platform getPlatform(){return this.platform;}
	public void setPlatform(final CPlatformType.Platform platform){this.platform=platform;}
	
	public DefaultExecutor createExecutor()
	{
		DefaultExecutor executor = new DefaultExecutor();
		if (this.workingDir!=null)
			executor.setWorkingDirectory(this.workingDir);
		if (this.exitValue!=null)
			executor.setExitValue(this.exitValue);
		if (this.watchdog!=null)
			executor.setWatchdog(new ExecuteWatchdog(this.watchdog));
		return executor;
	} 
	
	@Override
	public String toString()
	{
		return this.commandLine.toString();
	}
	
	public static class Output
	{
		protected int exitValue;
		protected String out;
		protected String err;
				
		public int getExitValue(){return this.exitValue;}
		public void setExitValue(final int exitValue){this.exitValue=exitValue;}
		
		public String getOut(){return this.out;}
		public void setOut(final String out){this.out=out;}

		public String getErr(){return this.err;}
		public void setErr(final String err){this.err=err;}
		
		public boolean hasOutput()
		{
			return CStringHelper.hasContent(this.out);
		}
		
		public boolean hasError()
		{
			return CStringHelper.hasContent(this.err);
		}
		
		public void dump()
		{
			System.out.println("exitValue: "+getExitValue());
			System.out.println("stdout: ["+getOut()+"]");
			System.out.println("stderr: ["+getErr()+"]");
		}
	}
}
