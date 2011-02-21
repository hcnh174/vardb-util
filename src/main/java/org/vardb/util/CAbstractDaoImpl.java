package org.vardb.util;

import java.io.Serializable;
import java.math.BigInteger;
import java.sql.Connection;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import javax.sql.DataSource;

import org.hibernate.Query;
import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.jdbc.Work;
import org.hibernate.stat.Statistics;
import org.springframework.beans.factory.annotation.Required;
import org.springframework.jdbc.core.JdbcTemplate;
import org.springframework.jdbc.datasource.SingleConnectionDataSource;
import org.springframework.stereotype.Repository;

@Repository
public abstract class CAbstractDaoImpl //implements IAbstractDao
{
	public static final int BATCH_SIZE=50;
	protected SessionFactory sessionFactory;
	private JdbcTemplate jdbcTemplate=null;
	
	@Required public void setSessionFactory(final SessionFactory sessionFactory){this.sessionFactory=sessionFactory;}
	
	public Session getSession()
	{
		return this.sessionFactory.getCurrentSession();
	}
	
	public Object load(Class<?> cls, Serializable id)
	{
		return getSession().load(cls,id);
	}
	
	public Object get(Class<?> cls, Serializable id)
	{
		return getSession().get(cls,id);
	}
	
	public void save(Object obj)
	{
		getSession().save(obj);
	}
		
	public void saveOrUpdate(Object obj)
	{
		getSession().saveOrUpdate(obj);
	}
	
	public void update(Object obj)
	{
		getSession().update(obj);
	}
	
	public void merge(Object obj)
	{
		getSession().merge(obj);
	}
	
	public void delete(Object obj)
	{
		getSession().delete(obj);
	}
	
	//////////////////////////////////////
	
	@SuppressWarnings("unchecked")
	public List find(StringBuilder buffer)
	{
		return find(buffer.toString());
	}
	
	@SuppressWarnings("unchecked")
	public List find(String sql)
	{
		System.out.println("FIND SQL="+sql);
		Query query=getSession().createQuery(sql);
		return query.list();
	}
	
	////////
	
	@SuppressWarnings("unchecked")
	public List findByNamedParam(StringBuilder buffer, Object... args)
	{
		return findByNamedParam(buffer.toString(),args);
	}
	
	@SuppressWarnings("unchecked")
	public List findByNamedParam(String sql, Object... args)
	{
		Map<String,Object> map=CStringHelper.createMap(args);
		Query query=getSession().createQuery(sql);
		for (String name : map.keySet())
		{
			Object value=map.get(name);
			//System.out.println("setting parameter: "+name+"="+value);
			query.setParameter(name,value);
		}
		//query.setProperties(map);
		return query.list();
	}
	
	////////
	
	public Object findUnique(StringBuilder buffer, Object... args)
	{
		return findUnique(buffer.toString(),args);
	}
	
	@SuppressWarnings("unchecked")
	public Object findUnique(String sql, Object... args)
	{
		List list=findByNamedParam(sql,args);
		if (list.size()>1)
			throw new CException("findUnique returned more than 1 result");
		if (list.size()==1)
			return list.get(0);
		else return null;
	}
	
	////////
	
	public void updateAll(Collection<? extends Object> list)
	{
		Session session=getSession();
		for (Object obj : list)
		{
			session.update(obj);
		}
		session.flush();
		session.clear();
	}
	
	public void saveOrUpdateAll(Collection<? extends Object> list)
	{
		if (list.isEmpty())
		{
			System.out.println("saveOrUpdateAll: list is empty");
			return;
		}
		Session session=getSession();
		for (Object obj : list)
		{
			session.saveOrUpdate(obj);
		}
		session.flush();
		session.clear();
	}
	
	public void saveAll(Collection<? extends Object> list)
	{
		if (list.isEmpty())
		{
			System.out.println("saveAll: list is empty");
			return;
		}
		Session session=getSession();
		for (Object obj : list)
		{
			session.save(obj);
		}
		session.flush();
		session.clear();
	}
	
	public void mergeAll(Collection<? extends Object> list)
	{
		if (list.isEmpty())
			return;
		Session session=getSession();
		for (Object obj : list)
		{
			session.merge(obj);
		}
		session.flush();
		session.clear();
	}
	
	public void deleteAll(Collection<? extends Object> list)
	{
		if (list.isEmpty())
			return;
		for (Object obj : list)
		{
			getSession().delete(obj);
		}
	}
	
	public void flush()
	{
		getSession().flush();
	}
	
	public void clear()
	{
		getSession().clear();
	}
	
	public void vacuuum()
	{
		String sql="VACUUM ANALYZE;\n";
		System.out.println("sql="+sql);
		getJdbcTemplate().execute(sql);
	}
	
	public void setDataSource(DataSource dataSource)
	{
		jdbcTemplate = new JdbcTemplate(dataSource);
    }
	
	public JdbcTemplate getJdbcTemplate()
	{
		if (jdbcTemplate==null)
			jdbcTemplate=createJdbcTemplate();
		return jdbcTemplate;
	}
	
	protected void getTotal(String sql, IPaging paging)
	{
		if (!paging.getPaged())
			return;
		Query query=getSession().createSQLQuery(sql);
		BigInteger total=(BigInteger)(query.uniqueResult());
		paging.setTotal(Integer.valueOf(total.toString()));
		System.out.println("total="+paging.getTotal());
	}
	
	private JdbcTemplate createJdbcTemplate()
	{
		getSession().doWork(new Work()
		{
			public void execute(Connection connection)
			{
				SingleConnectionDataSource dataSource=new SingleConnectionDataSource(connection,true);
				setDataSource(dataSource);
			}
		});
		return jdbcTemplate;
	}
	
	protected void getStatistics()
	{
		Statistics stats=sessionFactory.getStatistics();
		stats.logSummary();
		System.out.println("entity insert count: "+stats.getEntityInsertCount());
		System.out.println("entity load count: "+stats.getEntityLoadCount());
		System.out.println("entity update count: "+stats.getEntityUpdateCount());
		System.out.println("entity delete count: "+stats.getEntityDeleteCount());
		System.out.println("toString: "+stats.toString());
	}
}
