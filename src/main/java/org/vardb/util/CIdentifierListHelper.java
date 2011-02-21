package org.vardb.util;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public final class CIdentifierListHelper
{
	private static final int MAX_LIST_SIZE=100;	
	
	////////////////////////////////////////////////////////////////////////////

	private CIdentifierListHelper(){}
	
	public static String toSqlEnum(String field, Collection<? extends Enum<?>> list)
	{
		return toSql(field,CStringHelper.getNames(list));
	}
	
	public static String toSql(String field, Collection<String> list)
	{
		if (list.isEmpty())
			return "1=2";
		if (list.size()<MAX_LIST_SIZE)
			return createInList(field,list);
			//return CStringHelper.parenthesize(createInList(field,list));
			//return "("+createInList(field,list)+")";
		List<String> subqueries=new ArrayList<String>();
		createInListSubquery(field,list,subqueries);
		String sql=CStringHelper.join(subqueries," or ").trim();
		sql=CStringHelper.parenthesize(sql);
		return sql;
	}
	
	private static Collection<Collection<String>> split(Collection<String> ids)
	{
		Collection<Collection<String>> lists=new ArrayList<Collection<String>>();
		if (ids.size()<=MAX_LIST_SIZE)
		{
			lists.add(ids);
			return lists;
		}
		Collection<String> list=new ArrayList<String>();
		lists.add(list);
		for (String id : ids)
		{
			if (list.size()>=MAX_LIST_SIZE)
			{
				list=new ArrayList<String>();
				lists.add(list);
			}
			list.add(id);
		}
		return lists;
	}
	
	private static void createInListSubquery(String field, Collection<String> ids, List<String> subqueries)
	{
		// split the ids into smaller groups
		Collection<Collection<String>> lists=split(ids);
		for (Collection<String> list : lists)
		{
			addSubquery(subqueries,createInList(field,list));
		}
	}
	
	private static String createInList(String field, Collection<String> ids)
	{
		if (ids.isEmpty())
			return "";
		if (ids.size()==1)
			return field+"='"+ids.iterator().next()+"'";
		else return field+" in ('"+CStringHelper.join(ids,"','")+"')";
	}
	
	private static void addSubquery(List<String> subqueries, String subquery)
	{
		if (CStringHelper.hasContent(subquery))
			subqueries.add(subquery);
	}
}