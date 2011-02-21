package org.vardb.util;

import java.awt.image.BufferedImage;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Reader;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.imageio.ImageIO;
import javax.servlet.http.HttpServletRequest;

import org.apache.commons.httpclient.Header;
import org.apache.commons.httpclient.HttpClient;
import org.apache.commons.httpclient.HttpException;
import org.apache.commons.httpclient.HttpMethod;
import org.apache.commons.httpclient.HttpStatus;
import org.apache.commons.httpclient.NameValuePair;
import org.apache.commons.httpclient.cookie.CookiePolicy;
import org.apache.commons.httpclient.methods.GetMethod;
import org.apache.commons.httpclient.methods.PostMethod;
import org.apache.commons.httpclient.methods.multipart.FilePart;
import org.apache.commons.httpclient.methods.multipart.MultipartRequestEntity;
import org.apache.commons.httpclient.methods.multipart.Part;
import org.apache.commons.httpclient.methods.multipart.StringPart;
import org.apache.commons.io.IOUtils;
import org.vardb.util.CException;
import org.vardb.util.CFileHelper;
import org.vardb.util.CStringHelper;
import org.springframework.web.multipart.MultipartHttpServletRequest;
import org.springframework.web.multipart.commons.CommonsMultipartFile;

public class HttpHelper
{	
	private HttpHelper(){}
	
	public static String getRequest(String url, Map<String,Object> model)
	{
		CHttpRequest request=new CHttpRequest();
		return request.getRequest(url,model);
	}
	
	public static String getRequest(String url, Object... args)
	{
		Map<String,Object> params=CStringHelper.createMap(args);
		return getRequest(url,params);
	}
	
	public static String postRequest(String url, Map<String,Object> params)
	{
		CPostHttpRequest request=new CPostHttpRequest();
		return request.getRequest(url,params);
	}
	
	public static String postRequest(String url, Object... args)
	{
		Map<String,Object> params=CStringHelper.createMap(args);
		CPostHttpRequest request=new CPostHttpRequest();
		return request.getRequest(url,params);
	}
	
	public static BufferedImage postToImageGenerator(String url, Map<String,Object> model)
	{
		CPostHttpRequest request=new CPostHttpRequest();
		return request.postToImageGenerator(url,model);
	}	
	
	public static String postRedirectRequest(String url, Map<String,Object> model)
	{
		CPostHttpRequest request=new CPostHttpRequest();
		return request.getRequest(url,model);
	}
	
	public static String postMultipartRequest(String url, Map<String,Object> model)
	{
		Map<String,Object> files=Collections.emptyMap();//new HashMap<String,Object>();
		return postMultipartRequest(url,model,files);
	}
	
	public static String postMultipartRequest(String url, Map<String,Object> model, Map<String,Object> files)
	{
		PostMethod method=null;
		try
		{
			HttpClient client = new HttpClient();
			method = new PostMethod(url);        
			Part[] parts=new Part[model.size()+files.size()];
			int index=0;
			for (String name : model.keySet())
			{
				String value=(String)model.get(name);
				parts[index]=new StringPart(name,value);
				index++;
			}
			for (String name : files.keySet())
			{
				String filename=(String)files.get(name);
				File file=new File(filename);
				parts[index]=new FilePart(name,file);
				index++;
			}

			method.setRequestEntity(new MultipartRequestEntity(parts,method.getParams()));
			method.getParams().setCookiePolicy(CookiePolicy.BROWSER_COMPATIBILITY);

			// Execute the method.
			int statusCode = client.executeMethod(method);
			if (statusCode == HttpStatus.SC_OK)
			{
				return CHttpRequest.getResponseBody(method);
			}
			else if (statusCode == HttpStatus.SC_MOVED_TEMPORARILY)
			{
				Header header=method.getResponseHeader("Location");
				return header.getValue();
			}
			else
			{
				CStringHelper.println("Method failed: " + method.getStatusLine());
				return method.getStatusLine().toString();
			}
		}
		catch (HttpException e)
		{
			CStringHelper.println("Fatal protocol violation: " + e.getMessage());
			e.printStackTrace();
			throw new CException(e);
		}
		catch (IOException e)
		{
			CStringHelper.println("Fatal transport error: " + e.getMessage());
			e.printStackTrace();
			throw new CException(e);
		}
		finally
		{
			// Release the connection.
			method.releaseConnection();
		}
	}
		
	public static class CHttpRequest
	{
		protected Integer statusCode;
		
		public String getRequest(String url)
		{
			return getRequest(url,new HashMap<String,Object>());
		}

		public String getRequest(String url, Map<String,Object> model)
		{
			HttpClient client = new HttpClient();
			HttpMethod method = createMethod(url,model);
	        method.getParams().setCookiePolicy(CookiePolicy.BROWSER_COMPATIBILITY);
			try
			{
				// Execute the method.
				this.statusCode = client.executeMethod(method);
				if (this.statusCode==HttpStatus.SC_OK)
					return onOk(method);
				else if (this.statusCode==HttpStatus.SC_MOVED_TEMPORARILY)
					return onRedirect(method);
				else return onOther(method);
			}
			catch (HttpException e)
			{
				throw new CException(e);
			}
			catch (IOException e)
			{
				throw new CException(e);
			}
			finally
			{
				method.releaseConnection();
			}
		}
		
		public Integer getStatusCode(){return this.statusCode;}
		public void setStatusCode(Integer statusCode){this.statusCode=statusCode;}

		protected HttpMethod createMethod(String url, Map<String,Object> model)
		{
			// if there's already a query string, just add an &
			char separator=(url.indexOf('?')==-1) ? '?' : '&';
			url=url+separator+CWebHelper.createQueryString(model);
			return new GetMethod(url);
		}

		protected String onOk(HttpMethod method)
		{
			return getResponseBody(method);
		}
		
		protected String onRedirect(HttpMethod method)
		{
			CStringHelper.println("OnRedirect: " + method.getStatusLine());
			Header header=method.getResponseHeader("location");
			if (header==null)
				throw new CException("the response is invalid and did not provide the new location for the resource");
			return header.getValue();
		}
		
		protected String onOther(HttpMethod method)
		{
			String response=getResponseBody(method);
			CStringHelper.println("Method failed: " + method.getStatusLine()+", "+this.statusCode+", response="+response);
			return response;
		}
		
		protected static String getResponseBody(HttpMethod method)
		{
			try
			{
				InputStream input = method.getResponseBodyAsStream();
				return IOUtils.toString(input,CFileHelper.ENCODING.toString());
			}
			catch(IOException e)
			{
				throw new CException(e);
			}
		}
	}
	
	public static class CPostHttpRequest extends CHttpRequest
	{
		@Override
		protected HttpMethod createMethod(String url, Map<String,Object> model)
		{
			PostMethod method=new PostMethod(url);
			NameValuePair[] data = getNameValuePairs(model);
			method.setRequestBody(data);
			return method;
		}
		
		private static NameValuePair[] getNameValuePairs(Map<String,Object> model)
		{
			if (model==null)
				return new NameValuePair[0];
			NameValuePair[] data = new NameValuePair[model.size()];
			int index=0;
			for (String name : model.keySet())
			{
				//String value=(String)model.get(name);
				Object value=model.get(name);
				data[index]=new NameValuePair(name,value.toString());
				index++;
			}
			return data;
		}
		
		public BufferedImage postToImageGenerator(String url, Map<String,Object> model)
		{
			//Image img = Toolkit.getDefaultToolkit().createImage(resp.getData());			
			HttpClient client = new HttpClient();
			HttpMethod method = createMethod(url,model);
	        //method.getParams().setCookiePolicy(CookiePolicy.BROWSER_COMPATIBILITY);
			try
			{
				// Execute the method.
				this.statusCode = client.executeMethod(method);
				//String response=getResponseBody(method);
				if (this.statusCode==HttpStatus.SC_OK)
				{
					//byte[] data=method.getResponseBody();
					InputStream stream=method.getResponseBodyAsStream();
					return ImageIO.read(stream);
				}
				else throw new CException("returned statusCode: "+this.statusCode);
			}
			catch (HttpException e)
			{
				throw new CException(e);
			}
			catch (IOException e)
			{
				throw new CException(e);
			}
			finally
			{
				method.releaseConnection();
			}
		}
	}
	
	/////////////////////////////////////////////////////////////////////////////
	
	public static String simplePostRequest(String url)
	{
		return simplePostRequest(url,new HashMap<String,Object>());
	}
	
	public static String simplePostRequest(String path, Map<String,Object> model)
	{
		//System.out.println("posting request: path="+path);
		try
		{
			URL url=new URL(path);
			URLConnection con=url.openConnection();
			
			con.setDoOutput(true);
			OutputStream out=con.getOutputStream();
			OutputStream bout=new BufferedOutputStream(out);
			OutputStreamWriter writer=new OutputStreamWriter(bout);//,"8859_1");
			boolean first=true;
			for (String name : model.keySet())
			{
				String value=(String)model.get(name);
				//System.out.println("name="+name+", value="+value);
				if (!first)
				{
					writer.write("&");
					first=false;
				}
				writer.write(name+"="+value);
			}
			
			writer.flush();
			writer.close();
			
			InputStream stream=new BufferedInputStream(con.getInputStream());
			Reader reader=new BufferedReader(new InputStreamReader(stream));
			
			StringBuilder buffer=new StringBuilder();
			//int c;
			//while((c=reader.read())!=-1)
			for (int c=reader.read(); c!=-1; c=reader.read())
			{
				buffer.append((char)c);
			}
			return buffer.toString();
		}
		catch(MalformedURLException e)
		{
			throw new CException(e);
		}
		catch(IOException e)
		{
			throw new CException(e);
		}
	}
	
	////////////////////////////////////////////////////////////////////////////////
	
	public static void getFtpFiles(String server, String folder, List<String> filenames, String destination)
	{
		for (String filename : filenames)
		{
			String url=server+folder+filename;
			String outfile=destination+filename;
			getFtpFile(url,outfile);
		}
	}
	
	public static void getFtpFile(String url, String outfile)
	{
		BufferedInputStream in=null;
		BufferedOutputStream out=null;
		try
		{
			in = new BufferedInputStream(new URL(url).openStream());
			out = new BufferedOutputStream(new FileOutputStream(outfile),1024);
			byte[] data = new byte[1024];
			int x=0;
			while((x=in.read(data,0,1024))>=0)
			{
				out.write(data,0,x);
			}
		}
		catch(Exception e)
		{
			throw new CException(e);
		}
		finally
		{
			try
			{
				if (out!=null) out.close();
				if (in!=null) in.close();
			}
			catch (Exception e){}
		}
	}
	
	/*
	public static String getFtpFile(String server, String folder, String filename, String username, String password, String destination)
	{
		try
		{
			FTPClient ftp=openFTPConnection(server, username, password, folder);
			FTPFile[] files = ftp.listFiles();
			System.out.println("Number of files in dir: " + files.length);
			//DateFormat df = DateFormat.getDateInstance(DateFormat.SHORT);
			String file=null;
			for(int i=0; i<files.length; i++)
			{
				FTPFile ftpfile=files[i];
				if (!ftpfile.getName().equals(filename))
					continue;
				file=downloadFile(ftp,ftpfile,destination);
			}
			
			// Logout from the FTP Server and disconnect
			ftp.logout();
			ftp.disconnect();
			return file;
		}
		catch(Exception e)
		{
			throw new CException(e);
		}
	}
	
	
	public static List<String> getFtpFiles(String server, String folder, String regex, String username, String password, String destination)
	{
		try
		{
			//if (folder.charAt(folder.length()-1)!='/')
			//	folder+="/";
			FTPClient ftp=openFTPConnection(server, username, password, folder);
			FTPFile[] files = ftp.listFiles();
			Pattern pattern=Pattern.compile(regex);
			List<String> filenames=new ArrayList<String>();
			for(int i=0; i<files.length; i++)
			{
				FTPFile ftpfile=files[i];
				Matcher matcher=pattern.matcher(ftpfile.getName());
				if (!matcher.find())
					continue;
				filenames.add(downloadFile(ftp,ftpfile,destination));
			}
			// Logout from the FTP Server and disconnect
			ftp.logout();
			ftp.disconnect();
			return filenames;
		}
		catch(Exception e)
		{
			throw new CException(e);
		}
	}
	
	protected static FTPClient openFTPConnection(String server, String username, String password, String folder)
		throws IOException
	{
		FTPClient ftp=new FTPClient();
		ftp.connect(server);
		ftp.login(username, password);
		//System.out.println("Connected to " +  server + ".");
		//System.out.print(ftp.getReplyString());
		ftp.changeWorkingDirectory(folder);
		return ftp;
	}
	
	protected static String downloadFile(FTPClient ftp, FTPFile ftpfile, String destination)
		throws IOException
	{
		Date fileDate = ftpfile.getTimestamp().getTime();
		String destination_path=destination+File.separator+ftpfile.getName();
		File file = new File(destination_path);
		FileOutputStream fos = new FileOutputStream(file); 
		ftp.retrieveFile(file.getName(),fos);
		fos.close();
		if (file.setLastModified(fileDate.getTime()))
			CStringHelper.println("could not set last modified date on downloaded file");
		return destination_path;
	}
	*/
	
	/////////////////////////////////////////////////////////////////
	
	public static String getOriginalFilename(HttpServletRequest request, String name)
	{		
		if (!(request instanceof MultipartHttpServletRequest))
			throw new CException("request is not an instance of MultipartHttpServletRequest");	
		MultipartHttpServletRequest multipart=(MultipartHttpServletRequest)request;
		CommonsMultipartFile file=(CommonsMultipartFile)multipart.getFileMap().get(name);
		return file.getOriginalFilename();
	}
	
	/*
	public static String stripFiletype(String filename)
	{
		int index=filename.lastIndexOf('.');
		return filename.substring(0,index);
	}
	*/
}