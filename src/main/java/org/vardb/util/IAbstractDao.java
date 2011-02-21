package org.vardb.util;

import java.io.Serializable;
import java.util.Collection;
import java.util.List;

import org.hibernate.Session;
import org.springframework.jdbc.core.JdbcTemplate;

@SuppressWarnings("unchecked")
public interface IAbstractDao
{
	Object load(Class<?> cls, Serializable id);
	Object get(Class<?> cls, Serializable id);
	void save(Object obj);
	void saveOrUpdate(Object obj);
	void update(Object obj);
	void merge(Object obj);
	void delete(Object obj);	
	List find(StringBuilder buffer);
	List find(String sql);
	List findByNamedParam(StringBuilder buffer, Object... args);
	List findByNamedParam(String sql, Object... args);
	Object findUnique(StringBuilder buffer, Object... args);
	Object findUnique(String sql, Object... args);
	void updateAll(Collection<? extends Object> list);
	void saveOrUpdateAll(Collection<? extends Object> list);
	void saveAll(Collection<? extends Object> list);
	void mergeAll(Collection<? extends Object> list);
	void deleteAll(Collection<? extends Object> list);
	void flush();
	void clear();
	void vacuuum();
	Session getSession();
	JdbcTemplate getJdbcTemplate();
}
