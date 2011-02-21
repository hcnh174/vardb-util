package org.vardb.util.services;

public interface FreemarkerService
{
	String format(String path, Object... args);
	String formatStringTemplate(String str, Object... args);
}
