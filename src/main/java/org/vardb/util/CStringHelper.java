package org.vardb.util;

import java.io.PrintStream;
import java.io.UnsupportedEncodingException;
import java.net.URLDecoder;
import java.net.URLEncoder;
import java.nio.charset.Charset;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.text.DecimalFormat;
import java.text.Normalizer;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.commons.lang.builder.ReflectionToStringBuilder;
import org.apache.commons.lang.builder.ToStringBuilder;
import org.springframework.util.StringUtils;

import com.google.common.base.Charsets;
import com.google.common.base.Splitter;
import com.google.common.collect.Iterables;
import com.oreilly.servlet.Base64Decoder;
import com.oreilly.servlet.Base64Encoder;

public final class CStringHelper
{
	public static final String EMPTY_STRING="";
	public static final String UNICODE_SPACE=" ";
	public static final String ALPHABET="ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	
	private CStringHelper(){}
	
	public static String toString(Object obj)
	{
		return ReflectionToStringBuilder.reflectionToString(obj);
	}
	
	public static ToStringBuilder stringBuilder(Object obj)
	{
		return new ToStringBuilder(obj);
	}
	
	public static EqualsBuilder equalsBuilder()
	{
		return new EqualsBuilder();
	}
	
	public static HashCodeBuilder hashCodeBuilder()
	{
		return new HashCodeBuilder();
	}
	
	public static HashCodeBuilder hashCodeBuilder(int initialNonZeroOddNumber, int multiplierNonZeroOddNumber) 
	{
		return new HashCodeBuilder(initialNonZeroOddNumber,multiplierNonZeroOddNumber);
	}
	
	public static String trimAllWhitespace(String str)
	{
		return StringUtils.trimAllWhitespace(str);
	}
	
	/*
	public static String deleteWhitespace(String str)
	{
		return StringUtils.deleteWhitespace(str);
	}
	*/
	
	/*
	public static String stripWhitespace(String str)
	{
		//System.out.println("stripping whitespace: "+str);
		//return str.replace("^\\s+|\\s+$","");
		//return str.replace("\\s+","");
		StringBuilder buffer=new StringBuilder();
		for (int index=0;index<str.length();index++)
		{
			char ch=str.charAt(index);
			//if (Character.isLetterOrDigit(ch))
			if (!Character.isWhitespace(ch))
				buffer.append(ch);
		}
		return buffer.toString();
	}
	*/
	
	// removes HTML tags from string
	// probably better to use CRichTextFilter
	public static String stripHtml(String str)
	{
		return str.replaceAll("<\\/?[^>]+>","");
	}
	
	// trims each item in a collection
	// removes empty items
	public static List<String> trim(Collection<String> list)
	{
		List<String> ids=new ArrayList<String>();
		for (String id : list)
		{
			id=trim(id);
			if (!CStringHelper.isEmpty(id))
				ids.add(id);
		}
		return ids;
	}
	
	public static String replace(String str, String target, String replace)
	{
		return StringUtils.replace(str,target,replace);
		//Iterable<String> list=split(str,target);
		//return join(list,replace);
	}

	public static String join(Iterable<? extends Object> collection, String delimiter)
	{
		StringBuilder buffer=new StringBuilder();
		for (Iterator<? extends Object> iter=collection.iterator();iter.hasNext();)
		{
			buffer.append(iter.next().toString());
			if (iter.hasNext())
				buffer.append(delimiter);
		}
		return buffer.toString();
	}
		
	public static String join(Object[] array, String delimiter)
	{
		StringBuilder buffer=new StringBuilder();
		for (int index=0;index<array.length;index++)
		{
			buffer.append(array[index]);
			if (index<array.length-1)
				buffer.append(delimiter);
		}
		return buffer.toString();
	}
	
	public static String join(double[] array, String delimiter)
	{
		StringBuilder buffer=new StringBuilder();
		for (int index=0;index<array.length;index++)
		{
			buffer.append(array[index]);
			if (index<array.length-1)
				buffer.append(delimiter);
		}
		return buffer.toString();
	}

	// adds the "pad" character to the right as many times as necessary
	// to make the string the specified length
	// throws exception if truncated strings longer than the specified length
	public static String padRight(String str, char pad, int length)
	{
		int remainder=length-str.length();
		String padded=str+repeatString(String.valueOf(pad),remainder);
		if (padded.length()>length)
			throw new CException("padded string is longer than the specified length: ["+str+"] length="+length);
		return padded;
	}
	
	// adds the "pad" character to the left as many times as necessary
	// to make the string the specified length
	// throws exception if truncated strings longer than the specified length
	public static String padLeft(String str, char pad, int length)
	{
		int remainder=length-str.length();
		String padded=repeatString(String.valueOf(pad),remainder)+str;
		if (padded.length()>length)
			System.out.println("padded string is longer than the specified length: ["+str+"] length="+length);
			//throw new CException("padded string is longer than the specified length: ["+str+"] length="+length);
		return padded;
	}
	
	public static String repeatString(String str, int numtimes)
	{
		StringBuilder buffer=new StringBuilder("");
		for (int index=0;index<numtimes;index++)
		{
			buffer.append(str);
		}
		return buffer.toString();
	}
	
	public static String reverse(String str)
	{
		StringBuilder buffer=new StringBuilder(str);
		buffer.reverse();
		return buffer.toString();
	}
	
	public static String chunk(String original, int cols, String separator)
	{
		StringBuilder buffer=new StringBuilder();
		chunk(original,cols,separator,buffer);
		return buffer.toString();
	}
	
	public static void chunk(String original, int cols, String separator, StringBuilder buffer)
	{
		int length=original.length();
		int lines=(int)Math.ceil((double)length/(double)cols);
		int position=0;
		for (int index=0;index<lines;index++)
		{
			if (index<lines-1)
			{
				buffer.append(original.substring(position,position+cols));
				buffer.append(separator);
				position+=cols;
			}
			else buffer.append(original.substring(position));
		}
	}
	
	public static List<String> chunk(String original, int cols)
	{
		int length=original.length();
		int numlines=(int)Math.ceil((double)length/(double)cols);
		int position=0;
		List<String> lines=new ArrayList<String>();
		for (int index=0;index<numlines;index++)
		{
			if (index<numlines-1)
			{
				lines.add(original.substring(position,position+cols));
				position+=cols;
			}
			else lines.add(original.substring(position));
		}
		return lines;
	}
	
	public static Iterable<String> split(String raw, String delimiter, boolean clean)
	{
		if (!clean)
			return split(raw,delimiter);
		return clean(split(raw,delimiter));
	}
	
	public static Iterable<String> split(String raw, String delimiter)
	{
		if (raw==null)
			return Collections.emptyList();
		return Splitter.on(delimiter).split(raw);
	}
	
	public static List<String> splitAsList(String raw)
	{
		String delimiter=",";
		return splitAsList(raw,delimiter);
	}
	
	public static List<String> splitAsList(String raw, String delimiter)
	{
		List<String> list=new ArrayList<String>();
		Iterables.addAll(list,split(raw,delimiter));
		return list;
	}
	
	public static List<String> splitAsList(String raw, String delimiter, boolean clean)
	{
		List<String> list=new ArrayList<String>();
		Iterables.addAll(list,split(raw,delimiter));
		list=clean(list);
		return list;
	}
	
	public static String[] splitAsArray(String raw, String delimiter)
	{
		Iterable<String> iter=Splitter.on(delimiter).split(raw);
		return Iterables.toArray(iter, String.class);
	}
	
	/*
	public static List<String> split(String raw, String delimiter)
	{
		if (raw==null)
			return Collections.emptyList();
		String[] arr=raw.split(delimiter);
		List<String> list=new ArrayList<String>();
		list.addAll(Arrays.asList(arr));
		return list;
	}
	*/
	
	// splits on newlines, trims, and skips over blank lines
	public static List<String> splitLines(String str)
	{
		List<String> lines=new ArrayList<String>();
		for (String line : split(trim(str),"\n"))
		{
			line=trim(line);
			if (CStringHelper.isEmpty(line))
				continue;
			lines.add(line);
		}
		return lines;
	}
	
	public static List<Integer> splitInts(String str, String delimiter)
	{
		List<Integer> ints=new ArrayList<Integer>();
		if (!CStringHelper.hasContent(str))
			return ints;
		str=trim(str);		
		Iterable<String> list=split(str,delimiter);
		for (String item : list)
		{
			ints.add(Integer.valueOf(item));
		}
		return ints;
	}
	
	public static List<Float> splitFloats(String str, String delimiter)
	{
		List<Float> floats=new ArrayList<Float>();
		Iterable<String> list=split(str,delimiter);
		for (String item : list)
		{
			floats.add(Float.valueOf(item));
		}		
		return floats;
	}
	
	public static List<Double> splitDoubles(String str, String delimiter)
	{
		List<Double> doubles=new ArrayList<Double>();
		Iterable<String> list=split(str,delimiter);
		for (String item : list)
		{
			doubles.add(Double.valueOf(item));
		}		
		return doubles;
	}
	
	public static Collection<Collection<Integer>> split(Collection<Integer> ids, int max)
	{
		Collection<Collection<Integer>> lists=new ArrayList<Collection<Integer>>();
		Collection<Integer> list=new ArrayList<Integer>();
		lists.add(list);
		if (ids.size()<=max)
		{
			list.addAll(ids);
			return lists;
		}
		for (Integer id : ids)
		{
			if (list.size()>=max)
			{
				list=new ArrayList<Integer>();
				lists.add(list);
			}
			list.add(id);
		}
		return lists;
	}
	
	// trims each item and removes duplicates and empty lines
	public static List<String> clean(Iterable<String> items)
	{
		List<String> list=new ArrayList<String>();
		for (String item : items)
		{
			item=trim(item);
			if (!CStringHelper.isEmpty(item) && !list.contains(item))
				list.add(item);
		}
		return list;
	}
	
	public static List<String> clean(String[] items)
	{
		List<String> list=new ArrayList<String>();
		for (String item : items)
		{
			item=trim(item);
			if (!CStringHelper.isEmpty(item) && !list.contains(item))
				list.add(item);
		}
		return list;
	}
	
	public static List<String> trim(List<String> values)
	{
		List<String> trimmed=new ArrayList<String>();
		for (String value : values)
		{
			trimmed.add(trim(value));
		}
		return trimmed;
	}
	
	public static String[] trim(String[] values)
	{
		for (int index=0;index<values.length;index++)
		{
			values[index]=trim(values[index]);
		}
		return values;
	}
	
	public static String[] toArray(List<String> values)
	{
		String[] arr=new String[values.size()];
		values.toArray(arr);
		return arr;
	}
	
	
	public static String urlEncode(String str)
	{	
		return urlEncode(str,Charsets.US_ASCII);
	}
	
	public static String urlEncode(String str, Charset encoding)
	{
		try
		{
			return URLEncoder.encode(str,encoding.toString());
		}
		catch (UnsupportedEncodingException e)
		{
			throw new CException(e);
		}
	}
	
	public static String urlDecode(String str)
	{	
		return urlDecode(str,Charsets.US_ASCII);
	}
	
	public static String urlDecode(String str, Charset encoding)
	{
		try
		{
			return URLDecoder.decode(str,encoding.toString());
		}
		catch(UnsupportedEncodingException e)
		{
			throw new CException(e);
		}
	}
	
	public static String truncate(String str, int length)
	{
		if (str==null)
			return "";
		if (str.length()>=length)
			return str.substring(0,length);
		else return str;
	}
	
	public static String truncateEllipsis(String str, int length)
	{
		return truncate(str,length,"...");
	}
	
	// appends ellipsis or other trailing characters, adjusting the length accordingly
	public static String truncate(String str, int length, String trailing)
	{
		if (str==null)
			return "";
		if (str.length()<length)
			return str;
		int adjusted=length-trailing.length();
		return truncate(str,adjusted)+trailing;
	}
		
	public static String generateID()
	{
		try
		{
			// Initialize SecureRandom
			// This is a lengthy operation, to be done only upon
			// initialization of the application
			SecureRandom prng = SecureRandom.getInstance("SHA1PRNG");
	
			// generate a random number
			String randomNum = String.valueOf(prng.nextInt());
	
			// get its digest
			MessageDigest sha = MessageDigest.getInstance("SHA-1");
			byte[] result =  sha.digest(randomNum.getBytes());
			String id=hexEncode(result);
	
			//vSystem.out.println("Random number: " + randomNum);
			////System.out.println("Message digest: " + id);
			return id;
		}
		catch (NoSuchAlgorithmException e)
		{
			throw new CException(e);
		}
	}
	
	/**
	  * The byte[] returned by MessageDigest does not have a nice
	  * textual representation, so some form of encoding is usually performed.
	  *
	  * This implementation follows the example of David Flanagan's book
	  * "Java In A Nutshell", and converts a byte array into a String
	  * of hex characters.
	  *
	  * Another popular alternative is to use a "Base64" encoding.
	  */
	private static String hexEncode(byte[] input)
	{
		StringBuilder buffer=new StringBuilder();
		char[] digits = {'0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f'};
		for (int index = 0; index < input.length; ++index)
		{
			byte b = input[index];
			buffer.append(digits[(b&0xf0)>>4]);
			buffer.append(digits[b&0x0f]);
		}
		return buffer.toString();
	}
	
	public static String[] convertToArray(List<String> list)
	{
		String[] arr=new String[list.size()];
		list.toArray(arr);
		return arr;
	}
	
	public static String formatDecimal(Float value, int numdecimals)
	{
		if (value==null)
			return "";
		return formatDecimal((double)value,numdecimals);
	}
	
	public static String formatDecimal(Double value, int numdecimals)
	{
		if (value==null)
			return "";
		return String.format("%."+numdecimals+"f",value);
	}
	
	public static String formatScientificNotation(Double value, int numdecimals)
	{
		if (value==null)
			return "";
		String pattern="0."+CStringHelper.repeatString("#",numdecimals)+"E0";
		DecimalFormat format=new DecimalFormat(pattern);
		String formatted=format.format(value);
		if ("0E0".equals(formatted))
			formatted="0";
		return formatted;
	}
	
	public static void checkIdentifier(String identifier)
	{
		if (!CStringHelper.hasContent(identifier))
			throw new CException("identifier is null or empty: ["+identifier+"]");
	}
	
	public static boolean containsHtml(String str)
	{
		String regex="</?\\w+((\\s+\\w+(\\s*=\\s*(?:\".*?\"|'.*?'|[^'\">\\s]+))?)+\\s*|\\s*)/?>";
		return str.matches(regex);
	}
	
	public static boolean containsLinks(String str)
	{
		return str.toLowerCase().indexOf("<a href=")!=-1;
	}
	
	public static boolean isSpam(String str)
	{
		return (CStringHelper.containsHtml(str) || CStringHelper.containsLinks(str));
	}
	
	// slightly optimized to only trim if not null and length>0
	public static boolean hasContent(Object obj)
	{
		if (obj==null)
			return false;
		String value=obj.toString();
		if (value.length()==0)
			return false;
		value=trim(value);
		if (value.length()==0)
			return false;
		return !EMPTY_STRING.equals(value);
	}
	
	// skips null values
	@SuppressWarnings("unchecked")
	public static Map<String,Object> createMap(Object... args)
	{
		// if the only parameter was already a Map, return as is
		if (args.length==1 && args[0] instanceof Map)
			return (Map<String,Object>)args[0];
		int size=args.length/2;
		if (args.length%2!=0)
			throw new CException("name/value args should be provided as a multiple of 2: "+CStringHelper.join(args,","));
		Map<String,Object> map=new LinkedHashMap<String,Object>();
		for (int index=0;index<size;index++)
		{
			Object name=args[index*2];
			if (!(name instanceof String))
				throw new CException("parameter name at position "+index*2+" should be a String: "+CStringHelper.join(args,","));	
			Object value=args[index*2+1];
			if (name==null)
			{
				System.out.println("arg name is null for arg "+index+" and value "+value);
				continue;
			}
			if (value==null)
			{
				System.out.println("arg value is null for name "+name);
				continue;
			}
			map.put(name.toString(),value);
		}
		return map;
	}
	
	public static String extractBetween(String str, String prefix, String suffix)
	{
		int start=str.indexOf(prefix);
		if (start==-1)
			throw new CException("can't find prefix \""+prefix+"\" in string: "+str);
		start+=prefix.length();
		int end=str.indexOf(suffix,start);
		if (end==-1)
			throw new CException("can't find suffix \""+suffix+"\" in string: "+str);
		return str.substring(start,end);
	}

	public static String getText(byte[] bytes)
	{
		if (bytes==null)
			return null;
		StringBuilder buffer=new StringBuilder();
		for (int index=0;index<bytes.length;index++)
		{
			buffer.append((char)bytes[index]);
		}
		return buffer.toString();
	}
	
	public static String escapeSql(String value)
	{
		return CStringHelper.replace(value,"'","''");
	}
	
	public static List<String> escapeSql(Collection<String> values)
	{
		List<String> newvalues=new ArrayList<String>();
		for (String value : values)
		{
			newvalues.add(escapeSql(value));
		}
		return newvalues;
	}
	
	public static List<String> escapeSql(String[] values)
	{
		List<String> newvalues=new ArrayList<String>();
		for (String value : values)
		{
			newvalues.add(escapeSql(value));
		}
		return newvalues;
	}
	
	public static boolean isEmptyJson(String str)
	{
		return (!CStringHelper.hasContent(str) || "null".equals(str) || "{}".equals(str));
	}
	
	public static String encodeBase64(String unencoded)
	{
		return Base64Encoder.encode(unencoded);
	}
	
	public static String decodeBase64(String encoded)
	{
		return Base64Decoder.decode(encoded);
	}
	
	public static boolean isEmailAddress(String email)
	{
		return !(!CStringHelper.hasContent(email) || email.indexOf('@')==-1 || email.indexOf('.')==-1);
	}
	
	public static String parenthesize(String value)
	{
		return "("+value+")";
	}
	
	public static String quote(String str)
	{
		return doubleQuote(str);
	}
	
	public static String doubleQuote(String str)
	{
		return "\""+str+"\"";
	}
	
	public static String singleQuote(String str)
	{
		return "'"+str+"'";
	}
	
	public static String unquote(String str)
	{
		if (str.charAt(0)=='"' || str.charAt(0)=='\'')
			return str.substring(1,str.length()-1);
		return str;
	}
	
	public static String sqlQuote(String str)
	{
		return singleQuote(escapeSql(str));
	}
	
	public static boolean isEmpty(String str)
	{
		return str==null || EMPTY_STRING.equals(str);
	}
	
	public static List<String> wrap(Iterable<String> iter, String token)
	{
		return wrap(iter,token,token);
	}
	
	public static List<String> wrap(Iterable<String> iter, String before, String after)
	{
		List<String> items=new ArrayList<String>();
		for (String item : iter)
		{
			items.add(before+item+after);
		}
		return items;
	}
			
	
	// same as System.out.println, but uses Unicode
	//http://www.velocityreviews.com/forums/t137667-changing-system-out-encoding.html
	public static void println(String str)
	{
		PrintStream out = null;
		try
		{
			out = new PrintStream(System.err, true, CFileHelper.ENCODING.toString());
			out.println(str);
		}
		catch (UnsupportedEncodingException e)
		{
			System.err.println(e);
		}
	}
	
	public static void logError(String str)
	{
		System.err.println(str);
		if (CPlatformType.find().isWindows())
		{
			Date date=new Date();
			str=CDateHelper.format(date)+"\t"+str;
			CFileHelper.appendFile("c:/temp/errors.txt",str);
		}
	}
	
	// uses Springs StringUtils class, which uses Character.isWhitespace and should remove Japanese spaces as well
	public static String trim(String str)
	{
		return StringUtils.trimWhitespace(str);
	}
	
	public static int numOccurrences(String str, String target)
	{
		int count=0;
		int start=0;
		while ((start=str.indexOf(target,start))!=-1)
		{
			start+=target.length();
			count++;
		}
		return count;
	}
	

	public static String fixWideChars(String value)
	{
		if (value==null)
			return null;
		value=trim(value);
		value=normalize(value);
		value=fixWideNumbers(value);
		value=fixWideLetters(value);
		value=value.replace("ã€€"," ");
		value=value.replaceAll("  "," ");
		value=value.replace("?","?");
		value=value.replace("ã€œ","~");
		value=value.replace("ã€�",",");
		value=value.replace("ï¼�","/");
		return value;
	}
	
	public static String fixWideLetters(String value)
	{
		if (value==null)
			return null;
		String letters1="ï¼¡ï¼¢ï¼£ï¼¤ï¼¥ï¼¦ï¼§ï¼¨ï¼©ï¼ªï¼«ï¼¬ï¼­ï¼®ï¼¯ï¼°ï¼±ï¼²ï¼³ï¼´ï¼µï¼¶ï¼·ï¼¸ï¼¹ï¼º";
		String letters2="ABCDEFGHIJLKMNOPQRSTUVWXYZ";
		for (int index=0;index<letters1.length();index++)
		{
			String letter1=letters1.substring(index,index+1);
			String letter2=letters2.substring(index,index+1);
			//System.out.println("replacing "+letter1+" with "+letter2);
			value=value.replace(letter1,letter2);
		}
		return value;
	}
	
	public static String fixWideNumbers(String value)
	{
		if (value==null)
			return null;
		value=value.replace("ï¼‘","1");
		value=value.replace("ï¼’","2");
		value=value.replace("ï¼“","3");
		value=value.replace("ï¼”","4");
		value=value.replace("ï¼•","5");
		value=value.replace("ï¼–","6");
		value=value.replace("ï¼—","7");
		value=value.replace("ï¼˜","8");
		value=value.replace("ï¼™","9");
		value=value.replace("ï¼�","0");
		value=value.replace("ï¼Ž",".");
		return value;
	}
	
	public static String normalize(String value)
	{
		if (value==null)
			return null;
		value=Normalizer.normalize(value,Normalizer.Form.NFKC);
		//value=Normalizer.normalize(value,Normalizer.Form.NFD);
		return value;
	}
	
	public static List<String> getNames(Collection<? extends Enum<?>> items)
	{
		List<String> names=new ArrayList<String>();
		for (Enum<?> item : items)
		{
			names.add(item.name());
		}
		return names;
	}	
	
	/////////////////////////////////////
	
	public static int dflt(Integer value)
	{
		return dflt(value,0);
	}
	
	public static int dflt(Integer value, int dflt)
	{
		return (value==null) ? dflt : value;
	}	
	
	public static String dflt(String value)
	{
		return dflt(value,"");
	}
	
	public static String dflt(String value, String dflt)
	{
		return (value==null) ? dflt : value;
	}

	public static String dflt(Date value)
	{
		if (value==null)
			return "";
		return CDateHelper.format(value,CDateHelper.DATE_PATTERN);
	}
	
	public static String dflt(Object value)
	{
		return dflt(value,"");
	}
	
	public static String dflt(Object value, String dflt)
	{
		return (value==null) ? dflt : value.toString();
	}
	
	public static String joinNonEmpty(Iterable<String> iter, String delimiter)
	{
		List<String> items=clean(iter);
		return join(items,delimiter);
	}
	
	public static String getRandomWord(int min, int max)
	{
		int length=min+CMathHelper.randomInteger(max);
		StringBuilder buffer=new StringBuilder();
		for (int index=0;index<length;index++)
		{
			buffer.append(ALPHABET.charAt(CMathHelper.randomInteger(26)));
		}
		return buffer.toString();
	}
	
	/////////////////////////////////////
	
	// generates IDs
	public static void main(String[] argv)
	{
		int num=Integer.parseInt(argv[0]);
		for (int index=0;index<num;index++)
		{
			System.out.println(generateID());
		}
	}
}
