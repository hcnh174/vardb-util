package org.vardb.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.RandomAccessFile;
import java.io.Reader;
import java.io.StringReader;
import java.io.Writer;
import java.net.URL;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Properties;
import java.util.Scanner;
import java.util.zip.GZIPInputStream;

import org.xml.sax.InputSource;

import com.google.common.base.Charsets;
import com.google.common.io.Files;
import com.google.common.io.Resources;

public class CFileHelper
{
	public static final char SEP='/';
	public static final String SEPARATOR="/";
	public static final String NEWLINE="\n";
	public static final Charset ENCODING=Charsets.UTF_8;
	
	private CFileHelper(){}
	
	public static boolean exists(String filename)
	{
		File file=new File(filename);
		return file.exists();
	}

	public static Date getLastModifiedDate(String filename)
	{
		File file=new File(filename);
		return new Date(file.lastModified());
	}
	
	public static String readFile(String filename)
	{
		return readFile(filename,Charsets.UTF_8);
	}
	
	public static String readFile(String filename, Charset encoding)
	{
		checkPath(filename);
		File file=new File(filename);
		return readFile(file,encoding);
	}
	
	public static String readFile(File file, Charset encoding)
	{
		try
		{
			return Files.toString(file, encoding);	
		}
		catch (IOException e)
		{
			throw new CException(e);
		}
	}
	
	public static List<String> readLines(String filename, Charset encoding)
	{
		checkPath(filename);
		try
		{
			return Files.readLines(new File(filename), encoding);	
		}
		catch (IOException e)
		{
			throw new CException(e);
		}
	}
	
	public static List<String> listDirectories(String directoryname)
	{
		List<String> directories=new ArrayList<String>();

		try
		{
			File directory=new File(directoryname);
			if (directory.isFile())
				return null;
			String[] filenames=directory.list();

			for (int index=0;index<filenames.length;index++)
			{
				String filename=filenames[index];
				File file=new File(directoryname+filename);
				if (file.isDirectory())
					directories.add(filename);
			}
		}
		catch(Exception e)
		{
			throw new CException(e);
		}
		return directories;
	}

	public static List<String> listFiles(String dirname)
	{
		return listFiles(dirname,null);
	}
	
	public static List<String> listFiles(String dirname, String suffix)
	{
		return listFiles(dirname,suffix,false);
	}
	
	public static List<String> listFiles(String dirname, String suffix, boolean appendDir)
	{
		List<String> filenames=new ArrayList<String>();
		try
		{
			dirname=normalizeDirectory(dirname);
			File dir=new File(dirname);
			if (!dir.isDirectory())
				throw new CException("not a directory: ["+dirname+"]");
			File[] files=dir.listFiles();
			for (int index=0;index<files.length;index++)
			{
				File file=files[index];
				String filename=file.getName();
				if (suffix!=null && !filename.endsWith(suffix))
					continue;
				//System.out.println("filename="+filename);
				if (appendDir)
					filename=dirname+filename;
				filenames.add(filename);
			}
		}
		catch(Exception e)
		{
			throw new CException(dirname,e);
		}
		return filenames;
	}

	/*
	public static List<String> listFiles(String dirname, String suffix, boolean recursively)
	{
		if (recursively)
			return listFilesRecursively(dirname,suffix);
		else return listFiles(dirname,suffix);
	}
	*/
	
	public static List<String> listFilesRecursively(String dirname, String suffix, boolean recursively)
	{
		if (recursively)
			return listFilesRecursively(dirname,suffix);
		else return listFiles(dirname,suffix,true);
	}
	
	public static List<String> listFilesRecursively(String dirname)
	{
		List<String> filenames=new ArrayList<String>();
		listFilesRecursively(dirname,filenames,null,null);
		return filenames;
	}
	
	public static List<String> listFilesRecursively(String dirname, String suffix)
	{
		List<String> filenames=new ArrayList<String>();
		listFilesRecursively(dirname,filenames,suffix,null);
		return filenames;
	}
	
	public static List<String> listFilesRecursively(String dirname, String suffix, Date date)
	{
		List<String> filenames=new ArrayList<String>();
		listFilesRecursively(dirname,filenames,suffix,date);
		return filenames;
	}
	
	private static void listFilesRecursively(String dirname, List<String> filenames, String suffix, Date date)
	{
		try
		{
			dirname=normalizeDirectory(dirname);
			//System.out.println("dirname="+dirname);
			File dir=new File(dirname);
			if (!dir.isDirectory())
				throw new CException("not a directory");
			for (File file : dir.listFiles())
			{
				if (file.isDirectory())
					listFilesRecursively(file.getAbsolutePath(),filenames,suffix,date);
				String filename=file.getAbsolutePath();
				if (suffix!=null && !filename.endsWith(suffix))
					continue;
				//System.out.println("full filename="+filename+", lastModified="+file.lastModified()+" ("+new Date(file.lastModified())+")");
				if (date==null || file.lastModified()>=date.getTime())
				{
					//System.out.println("adding file filename="+filename+", lastModified="+file.lastModified()+" ("+new Date(file.lastModified())+")");
					//System.out.println("adding file filename="+filename+", lastModified="+CDateHelper.format(new Date(file.lastModified()))+")");
					filenames.add(filename);
				}
			}
		}
		catch(Exception e)
		{
			throw new CException(dirname,e);
		}
	}
	
	// returns false if any files could not be deleted
	public static boolean deleteFiles(String directoryname)
	{
		List<String> files=listFiles(directoryname);
		try
		{
			boolean result=true;
			for (int index=0;index<files.size();index++)
			{
				String filename=(String)files.get(index);
				File file=new File(filename);
				if (!file.delete())
				{
					CStringHelper.println("File could not be deleted: "+filename);
					result=false;
				}
			}
			return result;
		}
		catch(Exception e)
		{
			throw new CException(e);
		}
	}
	
	// returns true if the file was created
	public static boolean writeFileIfNotExists(String filename, String text)
	{
		if (exists(filename))
			return false;
		writeFile(filename,text);
		return true;
	}

	// creates an empty text file, useful for appending to
	public static void writeFile(String filename)
	{
		writeFile(filename,"");
	}
	
	public static void writeFile(String filename, String text)
	{
		writeFile(filename,text,ENCODING);
	}
	
	public static void writeFile(String filename, String text, Charset encoding)
	{
		try
		{
			checkPath(filename);
			Files.write(text, new File(filename), encoding);
		}
		catch(IOException e)
		{
			throw new CException(e);
		}
	}
	
	public static void appendFile(String filename, String str)
	{
		appendFile(filename,str,ENCODING,true);
	}
	
	public static void appendFile(String filename, String str, boolean newline)
	{
		appendFile(filename,str,ENCODING,newline);
	}
	
	public static void appendFile(String filename, String str, Charset encoding)
	{
		appendFile(filename,str,encoding,true);
	}
	
	public static void appendFile(String filename, String str, Charset encoding, boolean newline)
	{
		try
		{
			checkPath(filename);
			String end=newline ? "\n" : "";
			Files.append(str+end, new File(filename), encoding);
		}
		catch(IOException e)
		{
			throw new CException(e);
		}
	}

	public static void moveFile(String from, String to)
	{
		try
		{
			Files.move(new File(from), new File(to));
		}
		catch(Exception e)
		{
			throw new CException(e);
		}
	}
	
	public static void copyFile(String from, String to)
	{
		try
		{
			Files.copy(new File(from), new File(to));
		}
		catch(Exception e)
		{
			throw new CException(e);
		}
	}


	public static boolean deleteFile(String filename)
	{
		File file=new File(filename);
		
		boolean result=false;

		try
		{
			if (!file.isFile())
				return false;
			result=file.delete();
		}
		catch(Exception e)
		{
			throw new CException(e);
		}
		
		return result;
	}

	public static boolean deleteDirectory(String dirname)
	{
		File directory=new File(dirname);
		
		boolean result=false;

		try
		{
			if (!directory.isDirectory())
				return false;
			result=directory.delete();
		}
		catch(Exception e)
		{
			throw new CException(e);
		}
		
		return result;
	}
	
	/**
	Checks to see if the requested directory already exists - if not, creates it
	*/
	public static void createDirectory(String dir)
	{
		try
		{
			File path=new File(dir);
			if (path.exists())
				return;
			path.mkdir();
		}
		catch(SecurityException e)
		{
			throw new CException(e);
		}
	}

	public static String getResource(String path, Class<?> cls)
	{
		try
		{
			URL url=Resources.getResource(cls,path);
			return Resources.toString(url, ENCODING);
		}
		catch (IOException e)
		{
			throw new CException(e);
		}
	}
	
	// call without a leading slash, e.g. org/vardb/file.txt
	public static String getResource(String path)
	{
		try
		{
			URL url=Resources.getResource(path);
			return Resources.toString(url, ENCODING);
		}
		catch (IOException e)
		{
			throw new CException(e);
		}
	}
	
	// call without a leading slash, e.g. org/vardb/file.txt
	public static InputSource getResourceAsInputSource(String path)
	{
		InputStream stream=Thread.currentThread().getContextClassLoader().getResourceAsStream(path);
		InputStreamReader reader=new InputStreamReader(stream);
		return new InputSource(reader);
	}
	
	public static Properties getProperties(String path)
	{
		FileInputStream stream=null;
		try
		{
			Properties props=new Properties();
			stream=new FileInputStream(path);
			props.load(stream);
			return props;
		}
		catch (IOException e)
		{
			throw new CException(e);
		}
		finally
		{
			closeStream(stream);
		}
	}
	
	public static String getTempDirectory()
	{
		return normalizeDirectory(System.getProperty("java.io.tmpdir"));
	}
	
	public static String createTempFile(String prefix, String suffix)
	{
		try
		{
			File file=File.createTempFile(prefix,suffix);
			file.deleteOnExit();
			return file.getAbsolutePath();
		}
		catch(IOException e)
		{
			throw new CException(e);
		}
	}	
	
	public static String createTempFile(String prefix, String suffix, String text)
	{
		String filename=createTempFile(prefix,suffix);
		writeFile(filename,text);
		return filename;
	}
	
	public static String stripPath(String filepath)
	{
		filepath=normalize(filepath);
		int start=filepath.lastIndexOf(SEPARATOR);
		if (start==-1)
			start=0;
		else start+=1;
		return filepath.substring(start);
	}
	
	public static String stripFilename(String filepath)
	{
		filepath=normalize(filepath);
		int start=0;
		int end=filepath.lastIndexOf(SEPARATOR);
		if (end==-1)
			end=filepath.length();
		end+=1;
		return filepath.substring(start,end);
	}
	
	public static String stripFiletype(String filename)
	{
		int index=filename.lastIndexOf('.');
		return filename.substring(0,index);
	}
	
	public static int countLines(String filename)
	{
		Scanner scanner=createScanner(filename,"\n");
		int lines=0;
		while (scanner.hasNext())
		{
			scanner.next();
			lines++;
		}
		scanner.close();
		return lines;
	}
	
	public static int countLinesRecursively(String folder, String suffix, String outfile)
	{
		if (outfile!=null)
			CFileHelper.writeFile(outfile);
		int total_files=0;
		int total_lines=0;
		List<String> filenames=listFilesRecursively(folder,suffix);
		for (String filename : filenames)
		{
			total_files++;
			int count=countLines(filename);			
			total_lines+=count;
			System.out.println(filename+": "+count+" lines");
			if (outfile!=null)
				CFileHelper.appendFile(outfile,filename+"\t"+count);
		}
		System.out.println("total *"+suffix+" files: "+total_files);
		System.out.println("total *"+suffix+" lines: "+total_lines);
		return total_lines;
	}
	
	public static void main(String[] argv)
	{
		String folder=argv[0];
		String suffix=argv[1];
		String outfile=argv.length==3 ? argv[2] : null;
		CFileHelper.countLinesRecursively(folder,suffix,outfile);
	}
	
	public static String getDirFromFilename(String filename)
	{
		filename=normalize(filename);
		int end=filename.lastIndexOf(SEPARATOR);
		if (end==-1)
			end=0;
		else end+=1;
		return filename.substring(0,end);
	}

	public static String getIdentifierFromFilename(String filename, String suffix)
	{
		filename=normalize(filename);
		int start=filename.lastIndexOf(SEPARATOR);
		if (start==-1)
			start=0;
		else start+=1;
		int end=filename.indexOf(suffix);
		if (end==-1)
			throw new CException("can't find identifier in filename: "+suffix);
		return filename.substring(start,end);
	}
	
	public static String getIdentifierFromFilename(String filename)
	{
		return getIdentifierFromFilename(filename,getSuffix(filename));
	}
	
	public static String getSuffix(String filename)
	{
		filename=normalize(filename);
		int start=filename.lastIndexOf('.');
		if (start==-1)
			throw new CException("can't determine suffix from filename: "+filename);
		return filename.substring(start);
	}
	
	public static String changeSuffix(String filename, String oldsuffix, String newsuffix)
	{
		filename=normalize(filename);
		int start=filename.lastIndexOf(oldsuffix);
		if (start==-1)
			throw new CException("can't find suffix "+oldsuffix+" in filename: "+filename);
		return filename.substring(0,start)+newsuffix;
	}

	public static String getIdentifierFromDirname(String dir)
	{
		dir=normalize(dir);
		if (dir.charAt(dir.length()-1)==SEP)
			dir=dir.substring(0,dir.length()-1);
		int start=dir.lastIndexOf(SEPARATOR);
		if (start==-1)
			throw new CException("can't find path separator in dir: "+dir);
		start+=1;
		return dir.substring(start);
	}
	
	public static String normalize(String path)
	{
		path=path.replace("\\",SEPARATOR); // use back-slashes only
		path=path.replace("//",SEPARATOR); // remove double slashes￥
		path=path.replace("￥",SEPARATOR); // use back-slashes only
		return path;
	}
	
	public static String normalizeDirectory(String path)
	{
		String dir=normalize(path);
		if (dir.endsWith(SEPARATOR))
			return dir;
		return dir+SEPARATOR;
	}
	
	public static String getFullPath(String dir)
	{
		dir=normalizeDirectory(dir);
		File file=new File(dir);
		dir=file.getAbsolutePath();
		if (dir.endsWith("."));
			dir=dir.substring(0,dir.length()-1);
		return normalizeDirectory(dir);
	}
	
	private static void checkPath(String filename)
	{
		if ((filename.indexOf("c:")==0 || filename.indexOf("d:")==0) && !CPlatformType.find().isWindows())
				//CPlatformType.find().getPlatform()!=CPlatformType.Platform.WINDOWS)
			throw new CException("Windows path on a Unix platform: "+filename);
	}
	
	public static boolean isFolder(String path)
	{
		return normalize(path).endsWith(SEPARATOR);
	}
	
	public static String getTimestamp()
	{
		Date date=new Date();
		return String.valueOf(date.getTime());
	}
	
	public static void writeToOutputStream(OutputStream stream, String stdin)
	{
		try
		{
			if (stdin==null)
				return;
			System.out.println("writing stdin: "+stdin);
	        BufferedReader reader=new BufferedReader(new StringReader(stdin));
	        //String line;
	        //while ((line=reader.readLine())!=null)
	        for (String line=reader.readLine(); line!=null; line=reader.readLine())
			{
	        	stream.write(line.getBytes());
	        	stream.write("\n".getBytes());
			}
	        stream.flush();
	        stream.close();
		}
		catch (IOException t)
		{
			throw new CException(t);
		}
	}
	
	public static void closeReader(Reader reader)
	{
		if (reader==null)
			return;
		try
		{
			reader.close();
		}
		catch(IOException e)
		{
			System.out.println("error closing stream"+e.getMessage());
		}
	}
	
	public static void closeReader(RandomAccessFile reader)
	{
		if (reader==null)
			return;
		try
		{
			reader.close();
		}
		catch(IOException e)
		{
			System.out.println("error closing stream"+e.getMessage());
		}
	}
	
	public static void closeStream(InputStream stream)
	{
		if (stream==null)
			return;
		try
		{
			stream.close();
		}
		catch(IOException e)
		{
			System.out.println("error closing stream"+e.getMessage());
		}
	}
	
	public static void closeStream(OutputStream stream)
	{
		if (stream==null)
			return;
		try
		{
			stream.close();
		}
		catch(IOException e)
		{
			System.out.println("error closing stream"+e.getMessage());
		}
	}
	
	
	public static void closeWriter(Writer writer)
	{
		if (writer==null)
			return;
		try
		{
			writer.close();
		}
		catch(IOException e)
		{
			System.out.println("error closing writer"+e.getMessage());
		}
	}

	public static void unGzipFiles(String folder, List<String> filenames)
	{
		for (String filename : filenames)
		{
			unGzipFile(folder+filename);
		}
	}
	
	public static String unGzipFile(String infile)
	{
		GZIPInputStream in=null;
		OutputStream out=null;
		try
		{
	        // Open the compressed file
	        in = new GZIPInputStream(new FileInputStream(infile));
	        // Open the output file
	        String outfile=infile.substring(0,infile.lastIndexOf('.'));
	        out = new FileOutputStream(outfile);
	        // Transfer bytes from the compressed file to the output file
	        byte[] buf = new byte[1024];
	        int len;
	        while ((len = in.read(buf)) > 0)
	        {
	            out.write(buf, 0, len);
	        }
	        return outfile;
	    }
		catch (IOException e)
		{
			throw new CException("failed unGzipping file: "+infile,e);
	    }
		finally
		{
			// Close the file and stream
			 closeStream(in);
		     closeStream(out);
		}
	}
	
	public static Scanner createScanner(String filename)
	{
		return createScanner(filename,"\n");
	}
	
	public static Scanner createScanner(String filename, String delimiter)
	{
		try
		{
			FileReader reader=new FileReader(filename);
			return createScanner(reader,delimiter);
			//return new Scanner(reader).useDelimiter(delimiter);
		}
		catch(Exception e)
		{
			throw new CException(e);
		}
	}
	
	public static void closeScanner(Scanner scanner)
	{
		if (scanner==null)
			return;
		scanner.close();
	}
	
	public static Scanner createScanner(String filename, String delimiter, Charset encoding)
	{
		try
		{
			InputStreamReader reader=new InputStreamReader(new FileInputStream(filename),encoding.toString());
			return createScanner(reader,delimiter);
		}
		catch(Exception e)
		{
			throw new CException(e);
		}
	}
	
	public static Scanner createScanner(Reader in, String delimiter)
	{
		try
		{
			BufferedReader reader=new BufferedReader(in);
			return new Scanner(reader).useDelimiter(delimiter);
		}
		catch(Exception e)
		{
			throw new CException(e);
		}
	}
	
	public static String getFilenameAsUrl(String filename)
	{
		//return "file:"+filename;
		if (CPlatformType.find().isWindows())
			return "file:///"+filename;
		else return "file://"+filename;
	}
	
   public static String prompt(String message)
    {
		try
		{
			System.out.println(message);
			BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
			return br.readLine();
		}
		catch (IOException e)
		{
			throw new CException(e);
		}
	}
    
    public static boolean confirm(String message)
    {
    	String response=prompt(message);
    	return (response.equalsIgnoreCase("y") || response.equalsIgnoreCase("yes"));
    }
    
    public static String readInputStream(InputStream stream) throws IOException
    {
    	InputStreamReader reader=null;
    	try
    	{
			StringBuilder buffer=new StringBuilder();
			reader=new InputStreamReader(stream);
			for (int c=reader.read(); c!=-1; c=reader.read())
			{
				buffer.append((char)c);
			}
			return buffer.toString();
    	}
    	finally
    	{
    		CFileHelper.closeReader(reader);
    	}
    }
}