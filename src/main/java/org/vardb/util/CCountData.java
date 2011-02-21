package org.vardb.util;

import java.math.BigInteger;

public class CCountData
{	
	protected String name;
	protected String value;
	protected BigInteger count;
	
	public String getName(){return this.name;}
	public void setName(final String name){this.name=name;}

	public String getValue(){return this.value;}
	public void setValue(final String value){this.value=value;}

	public BigInteger getCount(){return this.count;}
	public void setCount(final BigInteger count){this.count=count;}
}