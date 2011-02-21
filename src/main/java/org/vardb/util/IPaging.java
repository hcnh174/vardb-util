package org.vardb.util;

public interface IPaging
{
	Integer getStart();
	void setStart(Integer start);

	Integer getPagesize();
	void setPagesize(Integer pagesize);

	int getTotal();
	void setTotal(int total);

	boolean getPaged();
	void setPaged(boolean paged);
	
	/*
	IFilter getFilter(); //hack!
	ISorting getSorting(); //hack!
	
	//hack!
	interface IFilter
	{
		String getFilter();
		void setFilter(final String filter);
	}
	*/
}
