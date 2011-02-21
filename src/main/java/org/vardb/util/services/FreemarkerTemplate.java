package org.vardb.util.services;

import java.util.Calendar;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.Id;
import javax.persistence.PreUpdate;
import javax.persistence.Table;
import javax.persistence.Transient;

@Entity
@Table(name="templates")
public class FreemarkerTemplate //implements IFreemarkerTemplate
{	
	private String name;
	private String content;
	private Calendar created;
	private Calendar updated;

	public FreemarkerTemplate(){}

	public FreemarkerTemplate(String name)
	{
		this.name=name;
	}
	
	@Id
	public String getName(){return this.name;}
	public void setName(String name){this.name=name;}

	public String getContent(){return this.content;}
	public void setContent(final String content){this.content=content;}
	
	@Column(insertable=false,updatable=false)
	public Calendar getCreated(){return this.created;}
	public void setCreated(final Calendar created){this.created=created;}

	@Column(insertable=false,updatable=false)
	public Calendar getUpdated(){return this.updated;}
	public void setUpdated(final Calendar updated){this.updated=updated;}
	
	@Transient
	public Calendar getLastModified(){return this.updated;}
	
	@PreUpdate
	protected void onUpdate()
	{
		this.updated = Calendar.getInstance();
		System.out.println("onUpdate called: "+this.updated.toString());
	}
}
