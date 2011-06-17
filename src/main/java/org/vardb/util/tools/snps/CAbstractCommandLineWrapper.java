package org.vardb.util.tools.snps;

import java.util.ArrayList;
import java.util.List;

import org.vardb.util.CFileHelper;

public abstract class CAbstractCommandLineWrapper
{
	protected IRuntimeService runtimeService=new RuntimeServiceImpl();
	protected String dir;
	protected String tempDir=null;
	protected List<String> tempfiles=new ArrayList<String>();
	protected String timestamp=null;
	protected boolean deleteTempfiles=true;
	
	public IRuntimeService getRuntimeService(){return this.runtimeService;}
	public void setRuntimeService(final IRuntimeService runtimeService){this.runtimeService=runtimeService;}
	
	public CAbstractCommandLineWrapper(String dir, String tempDir)
	{
		this.dir=dir;
		this.tempDir=tempDir;
	}
	
	public String createTempfile(String prefix, String suffix, String str)
	{
		String tempfile=this.tempDir+getRoot(prefix)+suffix;
		writeTempfile(tempfile,str);
		return tempfile;
	}
	
	public String expectOutfile(String prefix, String suffix)
	{
		String tempfile=this.tempDir+getRoot(prefix)+suffix;
		addTempfile(tempfile);
		return tempfile;
	}
	
	public String addTempfile(String tempfile)
	{
		this.tempfiles.add(tempfile);
		return tempfile;
	}
	
	public void writeTempfile(String tempfile, String str)
	{
		CFileHelper.writeFile(tempfile,str);
		addTempfile(tempfile);
	}
	
	public void deleteTempfiles()
	{
		if (!this.deleteTempfiles)
			return;
		for (String tempfile : this.tempfiles)
		{
			//System.out.println("deleting temp file "+tempfile);
			CFileHelper.deleteFile(tempfile);
		}
	}
	
	protected String getRoot(String prefix)
	{
		return prefix+"-"+getTimestamp();
	}
	
	protected String getTimestamp()
	{
		if (this.timestamp==null)
			this.timestamp=CFileHelper.getTimestamp();
		return this.timestamp;
	}
}