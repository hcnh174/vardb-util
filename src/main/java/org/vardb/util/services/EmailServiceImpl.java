package org.vardb.util.services;

import java.util.ArrayList;
import java.util.List;

import javax.annotation.Resource;

import org.springframework.beans.factory.annotation.Required;
import org.springframework.mail.MailSender;
import org.springframework.mail.SimpleMailMessage;
import org.springframework.stereotype.Service;
import org.springframework.transaction.annotation.Transactional;
import org.vardb.util.CException;
import org.vardb.util.CStringHelper;
import org.vardb.util.CWebHelper;

//@Service("emailService")
@Transactional(readOnly=false)
public class EmailServiceImpl implements EmailService
{
	@Resource(name="freemarkerService") private FreemarkerService freemarkerService;
	
	private MailSender mailSender;
	private List<String> emailHostnames=new ArrayList<String>();
	private String fromAddress;	
	
	public MailSender getMailSender(){return this.mailSender;}
	@Required public void setMailSender(MailSender mailSender){this.mailSender=mailSender;}
	
	public String getFromAddress(){return this.fromAddress;}
	@Required public void setFromAddress(final String fromAddress){this.fromAddress=fromAddress;}
	
	public void setEmailHostnames(String list)
	{
		for (String name : CStringHelper.split(list,",",true))
		{
			this.emailHostnames.add(name.toUpperCase());
		}
	}
	
	//////////////////////////////////////////////////
	
	public void sendEmail(String to, String subject, String body)
	{	
		SimpleMailMessage message=new SimpleMailMessage();
		message.setTo(to);
		message.setSubject(subject);
		message.setText(body);		
		sendEmail(message);
	}
	
	public void sendEmail(List<String> to, String subject, String body)
	{	
		SimpleMailMessage message=new SimpleMailMessage();
		message.setTo(CStringHelper.convertToArray(to));
		message.setSubject(subject);
		message.setText(body);		
		sendEmail(message);
	}
	
	public void sendEmail(String from, String to, String subject, String body)
	{	
		SimpleMailMessage message=new SimpleMailMessage();
		message.setFrom(from);
		message.setTo(to);
		message.setSubject(subject);
		message.setText(body);		
		sendEmail(message);
	}
	
	public void sendEmail(SimpleMailMessage message, String template, Object... args)//String name, Object value)
	{
		//String body=this.freemarkerService.format(template,name,value);
		String body=this.freemarkerService.format(template,args);
		System.out.println("template="+template);
		message.setText(body);
		sendEmail(message);
	}
	
	public void sendEmail(SimpleMailMessage message)
	{
		if (message.getTo().length==0)
			throw new CException("No To: addresses specified in email with subject: "+message.getSubject());
		try
		{
			if (!isEmailHostname())
				return;
			if (message.getFrom()==null)
				message.setFrom(this.fromAddress);
			this.mailSender.send(message);
		}
		catch(Exception e)
		{
			System.out.println(e);
			e.printStackTrace();
		}
	}

	private boolean isEmailHostname()
	{
		String server=CWebHelper.getServerName();
		return isEmailHostname(server);
	}
	
	private boolean isEmailHostname(String server)
	{
		return this.emailHostnames.contains(server.toUpperCase());
	}
}
