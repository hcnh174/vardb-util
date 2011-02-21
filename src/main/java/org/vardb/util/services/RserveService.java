package org.vardb.util.services;

import java.util.List;

import org.vardb.util.DataFrame;

public interface RserveService
{		
	String eval(final String str);	
	void eval(final List<String> str);
	void voidEval(final String str);
	List<Double> evalList(final String str);
	//DataFrame evalTable(final String str);
}