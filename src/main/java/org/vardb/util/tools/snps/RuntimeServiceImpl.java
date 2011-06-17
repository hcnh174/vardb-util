package org.vardb.util.tools.snps;

import java.io.ByteArrayOutputStream;
import java.io.IOException;

import org.apache.commons.exec.DefaultExecutor;
import org.apache.commons.exec.PumpStreamHandler;
import org.vardb.util.tools.snps.CCommandLine.Output;

public class RuntimeServiceImpl implements IRuntimeService
{
	public Output exec(CCommandLine commands)
	{	
		System.out.println("COMMAND LINE: "+commands.toString());
		DefaultExecutor executor=commands.createExecutor();
		Stream out = new Stream();
		Stream err = new Stream();
		executor.setStreamHandler(new PumpStreamHandler(out,err));
		Integer exitValue=null;
		try
		{
			exitValue=executor.execute(commands.getCommandLine());
		}
		catch (IOException e)
		{
			throw new CCommandLineException(e,commands.toString(),exitValue,out.getAsString(),err.getAsString());
		}
		CCommandLine.Output output=new CCommandLine.Output();
		output.setExitValue(exitValue);
		output.setOut(out.getAsString());
		output.setErr(err.getAsString());
		
		if (commands.getThrowExceptionIfOutputStream() && output.hasOutput())
			throw new CCommandLineException("Output stream is not empty",commands.toString(),output);
		if (commands.getThrowExceptionIfErrorStream() && output.hasError())
			throw new CCommandLineException("Error stream is not empty",commands.toString(),output);
		return output;
	}

	public static class Stream extends ByteArrayOutputStream
	{
		public String getAsString()
		{
			return new String(toByteArray()).trim();
		}
	}
}