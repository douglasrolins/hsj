package util;

public class SimException extends Exception
{

	/**
	 * 
	 */
	private static final long serialVersionUID = -6623181209901002495L;
	
	public SimException()
	{
		super ();
	}
	
	public SimException (String message)
	{
		super (message);
	}
	
	public SimException (String message, Throwable cause)
	{
		super (message, cause);
	}
	
	public SimException (String message, Throwable cause, boolean enableSuppression, boolean writableStackTrace)
	{
		super (message, cause, enableSuppression, writableStackTrace);
	}
	
	public SimException (Throwable cause)
	{
		super (cause);
	}
	
	
}
