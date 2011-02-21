package org.vardb.util.services;

import java.util.List;

import org.springframework.mail.SimpleMailMessage;

public interface EmailService
{	
	void sendEmail(String to, String subject, String body);
	void sendEmail(List<String> to, String subject, String body);
	void sendEmail(String from, String to, String subject, String body);
	void sendEmail(SimpleMailMessage message, String template, Object... args);
	void sendEmail(SimpleMailMessage message);	
}
