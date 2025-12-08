package util;

public class SimRuntimeException extends RuntimeException
{

	/**
	 * 
	 */
	private static final long serialVersionUID = -7116629342341012414L;
	
	public SimRuntimeException()
	{
		super ();
	}
	
	public SimRuntimeException (String message)
	{
		super (message);
	}
	
	public SimRuntimeException (String message, Throwable cause)
	{
		super (message, cause);
	}
	
	public SimRuntimeException (String message, Throwable cause, boolean enableSuppression, boolean writableStackTrace)
	{
		super (message, cause, enableSuppression, writableStackTrace);
	}
	
	public SimRuntimeException (Throwable cause)
	{
		super (cause);
	}
	

}
