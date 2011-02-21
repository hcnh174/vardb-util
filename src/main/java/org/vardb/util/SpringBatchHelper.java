package org.vardb.util;

import java.util.Map;

import org.springframework.batch.core.Job;
import org.springframework.batch.core.JobExecution;
import org.springframework.batch.core.JobParameters;
import org.springframework.batch.core.JobParametersBuilder;
import org.springframework.batch.core.JobParametersInvalidException;
import org.springframework.batch.core.launch.JobLauncher;
import org.springframework.batch.core.repository.JobExecutionAlreadyRunningException;
import org.springframework.batch.core.repository.JobInstanceAlreadyCompleteException;
import org.springframework.batch.core.repository.JobRestartException;
import org.springframework.context.ApplicationContext;
import org.springframework.context.support.ClassPathXmlApplicationContext;

public class SpringBatchHelper {
	
	public static JobExecution runJob(ApplicationContext context, String id, Object...args)
	{
		JobParameters jobparams=getJobParameters(args);
		
		JobLauncher jobLauncher = context.getBean(JobLauncher.class);
		Job job = (Job) context.getBean(id);
		
		try
		{
			JobExecution jobExecution = jobLauncher.run(job,jobparams);
			return jobExecution;
		}
		catch(JobParametersInvalidException e)
		{
			throw new CException(e);
		}
		catch(JobInstanceAlreadyCompleteException e)
		{
			throw new CException(e);
		}
		catch(JobExecutionAlreadyRunningException e)
		{
			throw new CException(e);
		}
		catch(JobRestartException e)
		{
			throw new CException(e);
		}
	}
	
	public static JobParameters getJobParameters(Object...args)
	{
		Map<String,Object> map=CStringHelper.createMap(args);
		JobParametersBuilder builder=new JobParametersBuilder();
		for (String name : map.keySet())
		{
			builder.addString(name,(String)map.get(name)); 
		}
		JobParameters jobparams=builder.toJobParameters();
		return jobparams;
	}
}
