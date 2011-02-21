package org.vardb.util;

public enum CPlatformType
{
	AIX("AIX",Platform.UNIX),
	DIGITAL_UNIX("Digital Unix",Platform.UNIX),
	FREE_BSD("FreeBSD",Platform.UNIX),
	HP_UX("HP UX",Platform.UNIX),
	IRIX("Irix",Platform.UNIX),
	LINUX("Linux",Platform.UNIX),
	MAC("Mac OS",Platform.UNIX),
	MPE("MPE/iX",Platform.OTHER),
	NETWARE("Netware 4.11",Platform.OTHER),
	OS2("OS/2",Platform.OTHER),
	SOLARIS("Solaris",Platform.UNIX),
	WIN2000("Windows 2000",Platform.WINDOWS),
	WIN95("Windows 95",Platform.WINDOWS),
	WIN98("Windows 98",Platform.WINDOWS),
	WINNT("Windows NT",Platform.WINDOWS),
	WINXP("Windows XP",Platform.WINDOWS),
	WINVISTA("Windows NT (unknown)",Platform.WINDOWS);
	
	public enum Platform{WINDOWS,UNIX,OTHER};
	protected String osname;
	protected Platform platform;
		
	CPlatformType(String osname, Platform platform)
	{
		this.osname=osname;
		this.platform=platform;
	}
	
	public String getOsname(){return this.osname;}
	public Platform getPlatform(){return this.platform;}
	
	public boolean isWindows(){return this.platform==Platform.WINDOWS;};
	public boolean isUnix(){return this.platform==Platform.UNIX;};
	
	public static CPlatformType find(String osname)
	{
		//System.out.println(osname);
		for (CPlatformType type : values())
		{
			if (type.getOsname().equals(osname))
				return type;
		}
		return null;
	}
	
	public static CPlatformType find()
	{
		String osname=System.getProperty("os.name");
		return find(osname);
	}
	
	public static String getBatchFileExtension()
	{
		CPlatformType.Platform platform=find().getPlatform();
		if (platform==CPlatformType.Platform.WINDOWS)
			return ".bat";
		else if (platform==CPlatformType.Platform.UNIX)
			return ".sh";
		else return ".sh";
	}
}
