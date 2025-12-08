package util;


//import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;

import org.apache.commons.math3.stat.descriptive.UnivariateStatistic;
import org.apache.commons.math3.stat.descriptive.moment.Mean;


public final class Simcalc
{
    public static final int VALUETYPE_INT       = 1;
	public static final int VALUETYPE_STRING    = 2;
	public static final int VALUETYPE_BYTES     = 3;
	//public static final int VALUETYPE_DEWEYID   = 4;
    
    //public static int i=0;
	/*public static MessageDigest md = null;					// for bloom filter construction

	static
	{
		try
		{
			md = MessageDigest.getInstance("MD5");
		}
		catch(NoSuchAlgorithmException ex)
		{
			return new XTCdiagMgr.print("Simcalc [init] Could not initialize MD5 digest. No algorithm found.");
		}
	}*/
	
	public final static byte getByte0 (int value)
    // return byte 0 from (byte[3]-byte[2]-byte[1]-byte[0]) of value
    {
        return (byte)(value & 255);
    }

    public final static byte getByte1 (int value)
    // return byte 1 from (byte[3]-byte[2]-byte[1]-byte[0]) of value
    {
        return (byte)((value >> 8) & 255);
    }

    public final static byte getByte2 (int value)
    // return byte 2 from (byte[3]-byte[2]-byte[1]-byte[0]) of value
    {
        return (byte)((value >> 16) & 255);
    }

    public final static byte getByte3 (int value)
    // return byte 3 from (byte[3]-byte[2]-byte[1]-byte[0]) of value
    {
        return (byte)((value >> 24) & 255);
    }

    public final static int getInt (byte byte0)
    // return int value of (byte0)
    {
        int ibyte0 = byte0;

        if (ibyte0 < 0) ibyte0 += 256;

        return ibyte0;
    }

    public final static int getInt (byte byte1, byte byte0)
    // return int value of (byte1-byte0)
    {
        int ibyte1 = byte1;
        int ibyte0 = byte0;

        if (ibyte1 < 0) ibyte1 += 256;
        if (ibyte0 < 0) ibyte0 += 256;

        return ( (ibyte1 << 8) | ibyte0 );
    }

    public final static int getInt (byte byte2, byte byte1, byte byte0)
    // return int value of (byte2-byte1-byte0)
    {
        int ibyte2 = byte2;
        int ibyte1 = byte1;
        int ibyte0 = byte0;

        if (ibyte2 < 0) ibyte2 += 256;
        if (ibyte1 < 0) ibyte1 += 256;
        if (ibyte0 < 0) ibyte0 += 256;

        return ( (ibyte2 << 16) | (ibyte1 << 8) | ibyte0 );
    }

    public final static int getInt (byte byte3, byte byte2, byte byte1, byte byte0)
    // return int value of (byte3-byte2-byte1-byte0)
    {
        int ibyte3 = byte3;
        int ibyte2 = byte2;
        int ibyte1 = byte1;
        int ibyte0 = byte0;

        if (ibyte3 < 0) ibyte3 += 256;
        if (ibyte2 < 0) ibyte2 += 256;
        if (ibyte1 < 0) ibyte1 += 256;
        if (ibyte0 < 0) ibyte0 += 256;

        return ( (ibyte3 << 24) | (ibyte2 << 16) | (ibyte1 << 8) | ibyte0 );
    }

	public final static byte[] getBytes(int length, int value)
	{
		if ((length < 0) || (length > 4))
			throw new SimRuntimeException(String.format("SimCalc (getBytes): Arrays of length %s are not supported", length));   
		
		byte[] result = new byte[length];
		
		for (int i = 1; i <= length; i++)
		{
			result[length - i] = (byte) (value & 255);
			value = (value >> 8);
		}
			
		return result;
	}
	
    public final static byte[] getBytes (int value)
    {
    	if (value < 0) throw new SimRuntimeException("Simcalc [getBytes]: Method supports only values greater than 0.");
    	
    	int valueLength = -1;
    	
    	if      (value < 256l) 			valueLength = 1;                   // 2 pow 8
    	else if (value < 65536l) 		valueLength = 2;                 // 2 pow 16
    	else if (value < 16777216l) 	valueLength = 3;              // 2 pow 24
    	else                        	valueLength = 4;            // 2 pow 32
    	
    	byte[] bytelist = new byte[valueLength];
    	
    	if (valueLength == 1)
    	{
    		bytelist[0] = getByte0 (value);
    	}
    	else if (valueLength == 2)
    	{
    		bytelist[0] = getByte1 (value);
    		bytelist[1] = getByte0 (value);
    	}
    	else if (valueLength == 3)
    	{
    		bytelist[0] = getByte2 (value);
    		bytelist[1] = getByte1 (value);
    		bytelist[2] = getByte0 (value);
    	}
    	else if (valueLength == 4)
    	{
    		bytelist[0] = getByte3 (value);
    		bytelist[1] = getByte2 (value);
    		bytelist[2] = getByte1 (value);
    		bytelist[3] = getByte0 (value);
    	}
    	
    	return bytelist;
    } // getBytes
    
    public final static int getInt (byte[] bytelist)
    {
    	if (bytelist.length == 1)
    	{
    		return getInt (bytelist[0]);
    	}
    	else if (bytelist.length == 2)
    	{
    		return getInt (bytelist[0],
    				       bytelist[1]);
    	}
    	else if (bytelist.length == 3)
    	{
    		return getInt (bytelist[0],
    				       bytelist[1],
    				       bytelist[2]);
    	}
    	else if (bytelist.length == 4)
    	{
    		return getInt (bytelist[0],
    				       bytelist[1],
	   				       bytelist[2],
    				       bytelist[3]);
    	}
    	
    	throw new SimRuntimeException ("Simcalc [getInt]: Invalid bytelist length " + bytelist.length);
    } // getInt
    
	/*public final static byte[] getBytes (String str)
	{
	    byte[] bytelist = null;
	
	    try
	    {
	        bytelist = str.getBytes (ConstantUtility.DEFAULT_CHARSET);
	        //if (!str.equals(Simcalc.getString(bytelist)))
	        //{
	        	//System.out.println(str);
	        //}
	        return bytelist;
	    }
	    catch (java.io.UnsupportedEncodingException e)
	    {
	        throw new SimRuntimeException (String.format("Unsupported encoding for String \" %s \""));
	    }
	
	} // getBytes from string
	*/

	/*public final static String getString (byte[] bytelist)
	{
		if (bytelist == null) return "<Bytelist is null>";
		
	    String data = null;
	
	    try
	    {
	        data = new String (bytelist, ConstantUtility.DEFAULT_CHARSET);
	
	        return data;
	    }
	    catch (java.io.UnsupportedEncodingException e)
	    {
	        throw new SimRuntimeException("SimCalc (getString): Unsupported encoding for bytelist.");
	    }
	
	} // getString
	*/    
    public final static int compare (byte[] value1, byte[] value2, int valueType)
	// returns -1 if value1<value2, 0 if value1==value2, 1 if value1>value2
	// comparison of deweyIDs takes NOT account of docID ! 
    {

    	if (valueType == VALUETYPE_INT)
    	{
    		int intValue1 = getInt (value1);
    		int intValue2 = getInt (value2);
    		
    		if (intValue1 < intValue2)  return -1;
    		if (intValue1 == intValue2) return 0;
    		
    		return 1; // longValue1 > longValue2
    	} // if valueType INT

    	else if ((valueType == VALUETYPE_BYTES) || (valueType == VALUETYPE_STRING))
    	{
    		int minlength = (value1.length < value2.length) ? value1.length : value2.length;
    		for(int i = 0; i < minlength; i++)
    		{
    			if (value1[i] < value2[i]) return -1; // value1 is less
    			if (value1[i] > value2[i]) return +1; // value2 is less
    		}
    		if(value1.length == value2.length)
    			return 0;
    		if(value1.length < value2.length)
    			return -1;
    		return 1;
    	}
    	
    	// old code by hauie
//    	else if ((valueType == VALUETYPE_BYTES) || (valueType == VALUETYPE_STRING))
//    	{
//    		if (value1.length < value2.length) return -1; // value1 is shorter and less
//    		if (value1.length > value2.length) return +1; // value2 is shorter and less
//    		
//    		// now value1.length == value2.length
//    		for (int i=0; i < value1.length; i++)
//    		{
//    			if (value1[i] < value2[i]) return -1; // value1 is less
//    			if (value1[i] > value2[i]) return +1; // value2 is less
//    		} // for i
//    		
//    		// now value1.length == value2.length and all bytes are equals -> values are equal
//    		return 0;
//    	} // if valueType STRING or BYTES
    	
    	/*else if (valueType == VALUETYPE_DEWEYID)
    	{
            // System.out.println("compare: " + ++i);
    		return XTCdeweyID.compare(value1, value2);
            // return new XTCdeweyID(-1, value1).compareTo(new XTCdeweyID(-1, value2));
    	}*/
    		
    	throw new SimRuntimeException ("Simcalc [compare]: Unsupported valueType " + valueType + ".");
    } // compare
    
    public final static int comparePrefix(byte[] value1, byte[] value2, int valueType)
    {
    	byte[] pValue1 = value1;
    	byte[] pValue2 = value2;

    	if (value1.length > value2.length)
    	{
    		pValue1 = new byte[value2.length];
    		System.arraycopy(value1, 0, pValue1, 0, value2.length);
    	}	
    	else if (value1.length < value2.length)
    	{
    		pValue2 = new byte[value1.length];
    		System.arraycopy(value2, 0, pValue2, 0, value1.length);
    	}
    	
    	return compare(pValue1, pValue2, valueType);
    }
    
    public final static boolean isPrefix(byte[] value1, byte[] value2)
    {
    	if(value2.length < value1.length)
    		return false;
    	
    	for(int i = 0; i < value1.length; i++)
    	{
    		if(value1[i] != value2[i])
    			return false;
    	}
    	return true;
    } // isPrefix
    public final static long getLong(byte[] b) 
    { 

    	int i; 
    	long val = 0; 
    	int bits = 64; 
    	long tval; 

    	for(i=8; i-- > 0;) 
    	{ 
    		bits -= 8; 
    		tval = b[i]; 
    		tval = tval < 0 ? 256 + tval : tval; 
    		val |= tval << bits; 
    	} 
    	return val; 
    } 
    public final static byte[] getBytes(long n)
    {
    	byte[] b = new byte[8];
    	b[7] = (byte) (n);
    	n >>>= 8;
    	b[6] = (byte) (n);
    	n >>>= 8;
    	b[5] = (byte) (n);
    	n >>>= 8;
    	b[4] = (byte) (n);
    	n >>>= 8;
    	b[3] = (byte) (n);
    	n >>>= 8;
        b[2] = (byte) (n);
        n >>>= 8;
        b[1] = (byte) (n);
        n >>>= 8;
        b[0] = (byte) (n);

        return b;
    }
    
	public final static long getLong(int upperInt, int lowerInt)
	{
		long value = ((lowerInt < 0) ? (long) (lowerInt + 4294967296L) : lowerInt);
		value = ((long) upperInt << 32) + value;
		return value; 
	} // getLong
	
	public final static int getUpperInt(long value)
	{
		return (int) (value >> 32);
	} // getUpperInt
	
	public final static int getLowerInt(long value)
	{
		return (int) value;
	} // getLowerInt
    
    /**
     * converts a boolean array to a byte array
     * @param bool boolean[]
     * @return byte[]
     * @author kschmidt
     **/
    public byte[] getBytes(boolean[] bool)
    {
      byte ret[] = new byte[bool.length / 8 + (bool.length % 8 >= 1 ? 1 : 0)];
      int arCnt = 0;
      for(int x = 0, y = 1; x != bool.length; x++, y*=2)
      {
        if(y > 128)
        {
          y = 1;
          arCnt++;
        }
        ret[arCnt] |= (byte) (y & (bool[x] ? y : 0));
      }
      return ret;
    }
    
    /**
     * returning a boolean array of a given byte array
     * @param values
     * @return boolean[]
     */

    public final static boolean[] getBoolean(byte[] values) {
    	boolean ret[] = new boolean[values.length*8];
    	int raise = 0;
    	for (int i=0; i< ret.length; i++)
    	{
    		raise = (int)Math.pow(2, i % 8);
    	      ret[i]= (values[i / 8] & raise) == raise;
    	}
    	return ret;
    }
    
    /**
     * Converts a boolean array to an integer. Note that if the array is greater than 32, only the last 32
     * bits will matter 
     *  
     */
    public final static int getInt(boolean[] a) //throws SimException
    {
        int n = 0, l = a.length;
        //if (l > 32)
        	//throw new SimException("Boolean array size cannot be greater than 32");
        for (int i = 0; i < l; ++i) 
        {
            n = (n << 1) + (a[i] ? 1 : 0);
        }
        return n;
    }

    
    
    
    /**
     * Sets bit <code>position</code> in <code>value</code> to <code>set</code>
     * @param value
     * @param position
     * @param set
     * @return
     */
    public final static byte setBit(int value, int position, boolean set)
    {
    	byte temp = 1;
    	temp = (byte) (temp << position);
    	
    	if (set)
    	{
    		return (byte) (temp | value);
    	}
    	else
    	{
    		return (byte) (~temp & value);
    	}
    } // setBit
    
    /**
     * Checks if bit <code>position</code> in <code>value</code> is set
     * @param value byte value to check
     * @param position postion of the bit that should be examined
     * @return
     */
    public final static boolean isBitSet(int value, int position)
    {
    	return (((value >>> position) & 1) > 0); 
    } // isSet
    
    /**
     * returns if boolean position in the byte array is true
     * @param value byte[]
     * @param boolIndex int
     * @return boolean
     * @author kschmidt
     **/
    public final static boolean isTrue(byte[] value, int boolIndex)
    {
      int raise = (int)Math.pow(2, boolIndex % 8);
      return (value[boolIndex / 8] & raise) == raise;
    }   
    
    
    //  Returns a bitset containing the values in bytes.
    // The byte-ordering of bytes must be big-endian which means the most significant bit is in element 0.
    public final static BitSet getBitSet(byte[] bytes) 
    {
        BitSet bits = new BitSet();
        for (int i=0; i<bytes.length*8; i++) {
            if ((bytes[bytes.length-i/8-1]&(1<<(i%8))) > 0) {
                bits.set(i);
            }
        }
        return bits;
    }

    public final static BitSet getBitSetCompensate(byte[] bytes) 
    {
        BitSet bits = new BitSet();
        for (int i=0; i<bytes.length*8; i++) {
            if ((bytes[bytes.length-i/8-1]&(1<<(i%8))) > 0) {
                bits.set(i);
            }
        }
        bits.set(bytes.length*8);
        return bits;
    }
    
    // Returns a byte array of at least length 1.
    // The most significant bit in the result is guaranteed not to be a 1
    // (since BitSet does not support sign extension).
    // The byte-ordering of the result is big-endian which means the most significant bit is in element 0.
    // The bit at index 0 of the bit set is assumed to be the least significant bit.
    public final static byte[] getBytes(BitSet bits) 
    {
        byte[] bytes = new byte[bits.length()/8+1];
        for (int i=0; i<bits.length(); i++) {
            if (bits.get(i)) {
                bytes[bytes.length-i/8-1] |= 1<<(i%8);
            }
        }
        return bytes;
    }

    public final static byte[] getBytesCompensate(BitSet bits) 
    {
        byte[] bytes = new byte[(bits.length()-1)/8+1];
        for (int i=0; i<(bits.length()-1); i++) {
            if (bits.get(i)) {
                bytes[bytes.length-i/8-1] |= 1<<(i%8);
            }
        }
        return bytes;
    }
    
    /*
     * Generates a bloom filter on the vocIds contained in the given path.
     */
	/*public final static int[] generateBloomFilter(Stack<Integer> path)
	{
		// int array to store the bloom filter values
		int[] bloomFilter = new int[XTCsvrCfg.fluxIndexBloomFilterSize];
		for(int i = 0; i < bloomFilter.length; i++)
			bloomFilter[i] = 0;
		
		byte[] bitNumberBytes = new byte[4];	// 4-byte partition of the message digest
		int    bitNumber  = 0;					// number of bit to set
		int    byteNumber = 0;					// number of byte to set
		byte[] message = null;					// md5 message
		for(Integer vocID : path)
		{
			// generate 128-bit MD5 message (16x8 bits)
			message = md.digest(Simcalc.getBytes(vocID));
		
			// partition MD5 message
			for(int i = 0; i < 4; i++) // setting four bits in the bloom filter
			//for(int i = 0; i < 2; i++) // setting one bit in the bloom filter
			{
				// create integer number from ith 4-bytes of the message
				System.arraycopy(message, i, bitNumberBytes, 0, bitNumberBytes.length);
				bitNumber = Simcalc.getInt(bitNumberBytes);
				if(bitNumber < 0)
					bitNumber = -bitNumber;
				
				// generates a positive bitNumber (0 .. bloomFilterSize-1)
				bitNumber %= XTCsvrCfg.fluxIndexBloomFilterSize*8;
				
				// set the required bit number in the bloom filter
				//System.out.println("turn : " + i + " setting bit position: " + bitNumber);
				byteNumber = (int)bitNumber/8;
				bitNumber  = bitNumber%8;
				//System.out.println("turn: " + i + " setting bit position: " + bitNumber + " in byte " + byteNumber);
				bloomFilter[byteNumber] |= (int)Math.pow(2, bitNumber);
			}
		}		
		return bloomFilter;
	}*/

	/**
	 * Matches two bloom filters.
	 * 
	 * @param qFilter the bloom filter generated for the query path
	 * @param iFilter the bloom filter from the indexed record
	 * @return
	 */
	public final static boolean matchBloomFilters(int[] qFilter, byte[] iFilter)
	{
		if(qFilter.length != iFilter.length)
			return false;
		
		int iFilterComponent = 0;
		for(int i = 0; i < qFilter.length; i++)
		{
			iFilterComponent = (int)iFilter[i];
			if((iFilterComponent | qFilter[i]) != iFilterComponent)
				return false;
		}
		return true;
	}
	
	public final static byte[] cropByteArray(byte[] array, int length)
	{
		if(length >= array.length)
			return array;
		
		byte[] result = new byte[length];
		for(int i = 0; i < length; i++)
			result[i] = array[i];
		return result;
	}
     
    // old code by cm 14.07.2006
//    public final static boolean isPrefix(byte[] value1, byte[] value2) throws XTCexception
//    {
//    	int minLength = (value1.length < value2.length) ? value1.length : value2.length;
//    	for(int i = 0; i < minLength; i++)
//    	{
//    		if(value1[i] != value2[i])
//    			return false;
//    	}
//    	return true;
//    } // isPrefix
	
	 //by lr on 26.03.07 begin of inclusion
  public static byte[] getBytes(double d)
  {
  	byte[] dBytes = new byte[8];
  	long l = Double.doubleToLongBits(d);
  	dBytes[0] = (byte)(l >>> 56);
  	dBytes[1] = (byte)(l >>> 48);
  	dBytes[2] = (byte)(l >>> 40);
  	dBytes[3] = (byte)(l >>> 32);
  	dBytes[4] = (byte)(l >>> 24);
  	dBytes[5] = (byte)(l >>> 16);
  	dBytes[6] = (byte)(l >>>  8);
  	dBytes[7] = (byte)(l >>>  0);
  	return dBytes;
  	
  	//byte[] bytes = new byte[8];
    //ByteBuffer.wrap(bytes).putDouble(d);
    //return bytes;
  }
  
  public static double getDouble(byte[] byteList)
  {
  	if (byteList.length != 8) throw new SimRuntimeException("Simcalc [getDouble]: Wrong byteList size");
  	long l = ((byteList[7] & 0xFFL) << 0) +
 	 	((byteList[6] & 0xFFL) << 8) +
 	 	((byteList[5] & 0xFFL) << 16) +
 	 	((byteList[4] & 0xFFL) << 24) +
 	 	((byteList[3] & 0xFFL) << 32) +
 	 	((byteList[2] & 0xFFL) << 40) +
 	 	((byteList[1] & 0xFFL) << 48) +
 	 	(((long) byteList[0]) << 56);
  	return Double.longBitsToDouble(l);
  	
  	//let's simplify the matters 
  	//return ByteBuffer.wrap(byteList).getDouble();
  }
  
  public static double getDouble(byte[] byteList, int offset)
  {
  	if (byteList.length < offset + 8) throw new SimRuntimeException("Simcalc [getDouble]: Invalid bytelist length=" + byteList.length);
  	//return ByteBuffer.wrap(byteList, offset, 8).getDouble();
  	long l = ((byteList[offset + 7] & 0xFFL) << 0) +
 	 	((byteList[offset + 6] & 0xFFL) << 8) +
 	 	((byteList[offset + 5] & 0xFFL) << 16) +
 	 	((byteList[offset + 4] & 0xFFL) << 24) +
 	 	((byteList[offset + 3] & 0xFFL) << 32) +
 	 	((byteList[offset + 2] & 0xFFL) << 40) +
 	 	((byteList[offset + 1] & 0xFFL) << 48) +
 	 	(((long) byteList[offset + 0]) << 56);
  	return Double.longBitsToDouble(l);
  }
  
  public static int getInt (byte[] byteList, int offset)
  {
  	if (byteList.length < offset + 4) throw new SimRuntimeException ("Simcalc [getInt]: Invalid bytelist length: " + byteList.length + " offset:" + offset + ")");
  	return getInt (byteList[offset + 0],
  			       byteList[offset + 1],
   			       byteList[offset + 2],
  			       byteList[offset + 3]);
  	
  }
  
	/*public static String getString (byte[] byteList, int offset, int size)
	{
		if (byteList.length < offset + size) throw new SimRuntimeException ("Simcalc [getInt]: Invalid bytelist length: " + byteList.length + " offset:" + offset + ")");
		byte[] byteString = new byte[size];
		System.arraycopy(byteList, offset, byteString, 0, size);
		return getString(byteString);
	}*/
  
	public final static byte[] putBytes(int value, byte[] byteList, int offset)
	{
		if (byteList.length < offset + 4) throw new SimRuntimeException ("No room for an integer in byte array (length:" + byteList.length + " offset:" + offset + ")");   
		
		byteList[offset + 3] = (byte) (value >>> 0);
		byteList[offset + 2] = (byte) (value >>> 8);
		byteList[offset + 1] = (byte) (value >>> 16);
		byteList[offset + 0] = (byte) (value >>> 24);	
		return byteList;
	}
	
	
	public final static byte[] putBytes(double value, byte[] byteList, int offset)
	{
		if (byteList.length < offset + 8) throw new SimRuntimeException ("No room for a double in byte array (length:" + byteList.length + " offset:" + offset + ")");   
		
		
		long l = Double.doubleToLongBits(value);
		byteList[offset + 0] = (byte)(l >>> 56);
		byteList[offset + 1] = (byte)(l >>> 48);
		byteList[offset + 2] = (byte)(l >>> 40);
		byteList[offset + 3] = (byte)(l >>> 32);
		byteList[offset + 4] = (byte)(l >>> 24);
		byteList[offset + 5] = (byte)(l >>> 16);
		byteList[offset + 6] = (byte)(l >>>  8);
		byteList[offset + 7] = (byte)(l >>>  0);
	  	
	    //ByteBuffer.wrap(byteList).putDouble(offset, value);
		return byteList;
	}
  
  public static short getShort (byte[] byteList, int offset)
  {
  	if (byteList.length < offset + 2) throw new SimRuntimeException ("Simcalc [getShort]: Invalid bytelist length=" + byteList.length);
  	return  (short) getInt (byteList[offset + 0],
  			                byteList[offset + 1]);
  	
  } //lr end of inclusion
  
  
  public static final List<Integer> getIntList (long l)
  {
	  final int chop = Integer.MAX_VALUE;
	  List<Integer> list = new ArrayList<Integer>();
	  while (l > chop)
	  {
		  list.add(chop);
		  l -= chop;
	  }
	  list.add((int)l);
	  return list;
  }
  
	public static double calcMean(double[] values, boolean removeExtremes)
	{
		UnivariateStatistic stat = new Mean();
		if (removeExtremes && values.length > 2)
		{
			Arrays.sort(values);
			return stat.evaluate(values, 1, values.length - 2);
		}
		else
		{
			return stat.evaluate(values);
		}
		
	}
	
	/*public static FractionRep decimalToFraction(double value, double accuracy)
	{
	    // ===================================================
	    // Richards implementation by Sjaak
	    //
		//https://stackoverflow.com/questions/5124743/algorithm-for-simplifying-decimal-to-fractions/42085412#42085412
		
	    // Split value in a sign, an integer part, a fractional part
	    int sign = value < 0 ? -1 : 1;
	    value = value < 0 ? -value : value;
	    int integerpart = (int)value;
	    value -= integerpart;
	
	    // check if the fractional part is near 0
	    double minimalvalue = value - accuracy;
	    if (minimalvalue < 0.0) return new FractionRep(sign * integerpart, 1);
	
	    // check if the fractional part is near 1
	    double maximumvalue = value + accuracy;
	    if (maximumvalue > 1.0) return new FractionRep(sign * (integerpart + 1), 1);
	
	    // Richards
	    double z = value;
	    int denominator0 = 0;
	    int denominator1 = 1;
	    int numerator0 = 1;
	    int numerator1 = 0;
	    int n = (int)z;
	    while (true)
	    {
	        z = 1.0 / (z - n);
	        n = (int)z;
	
	        int temp = denominator1;
	        denominator1 = denominator1 * n + denominator0;
	        denominator0 = temp;
	
	        temp = numerator1;
	        numerator1 = numerator1 * n + numerator0;
	        numerator0 = temp;
	
	        double d = (double)numerator1 / denominator1;
	        if (d > minimalvalue && d < maximumvalue) break;
	    }
	    return new FractionRep(sign * (integerpart * denominator1 + numerator1), denominator1);
	}
	*/
    
} // class Simcalc
