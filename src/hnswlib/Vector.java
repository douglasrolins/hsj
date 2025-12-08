package hnswlib;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import io.jhdf.HdfFile;
import io.jhdf.api.Dataset;
import io.jhdf.object.datatype.DataType;

public class Vector implements Item<Integer, double[]>
{

	private static final long serialVersionUID = 1L;

	private final int id;
	private final int foreingKey;
	private final double[] vector;

	public Vector(int id, int foreingKey, double[] vector)
	{
		this.id = id;
		this.foreingKey = foreingKey;
		this.vector = vector;
	}

	@Override
	public Integer id()
	{
		return id;
	}

	public int foreingKey()
	{
		return foreingKey;
	}

	@Override
	public double[] vector()
	{
		return vector;
	}

	@Override
	public int dimensions()
	{
		return vector.length;
	}

	@Override
	public String toString()
	{
		return "Vector{" + "id='" + id + '\'' + "foreignKey='" + foreingKey + '\'' + ", vector="
				+ Arrays.toString(vector) + '}';
	}

	public void normalize()
	{
		double norm = 0.0;
		for (double d : vector)
		{
			norm += d * d;
		}
		norm = Math.sqrt(norm);

		for (int i = 0; i < vector.length; i++)
		{
			vector[i] /= norm;
		}
	}

	public static void normalizeVectors(List<Vector> vectors)
	{
		for (Vector v : vectors)
		{
			v.normalize();
		}
	}

	public static List<Vector> fromCsvFile(Path sourcePath)
	{
		List<Vector> vectors = new ArrayList<>();

		try (BufferedReader reader = new BufferedReader(
				new InputStreamReader(new FileInputStream(sourcePath.toString()))))
		{

			int count = 1;
			String line = "";

			while ((line = reader.readLine()) != null)
			{
				String[] tokens = line.split(",");
				//String id = String.valueOf(count);
				int id = Integer.valueOf(count);

				//String foreingKey = tokens[0];
				int foreingKey = Integer.valueOf(tokens[0]);

				double[] vector = new double[tokens.length - 1];
				for (int i = 1; i < tokens.length - 1; i++)
				{
					vector[i] = Float.parseFloat(tokens[i]);
				}

				Vector w = new Vector(id, foreingKey, vector);
				vectors.add(w);
				count++;
			}
			//reader.close();
		} catch (Exception e)
		{
			System.out.print("Error occurred while reading csv file: " + e);
		}

		return vectors;
	}

	public static List<Vector> fromHdf5File(Path sourcePath, String datasetHdf5)
	{
		List<Vector> vectors = new ArrayList<>();
		try (HdfFile hdfFile = new HdfFile(Paths.get(sourcePath.toString())))
		{

			// set collection in hdf5 file
			Dataset dataset = hdfFile.getDatasetByPath(datasetHdf5);

			// get data type
			DataType datatype = dataset.getDataType();
			Class<?> dataTypeClass = datatype.getJavaType();
			
			// float
			if (dataTypeClass == float.class)
			{

				float[][] floatArray = (float[][]) dataset.getData();

				for (int i = 0; i < floatArray.length; i++)
				{

					// convert to double
					double[] doubleArray = new double[floatArray[i].length];
					for (int j = 0; j < floatArray[i].length; j++)
					{
						doubleArray[j] = (double) floatArray[i][j];
					}

					// new vector: id, foreingKey, vector
					vectors.add(new Vector(i, 0, doubleArray));
				}
				
				floatArray = null;

			}

			// double
			else if (dataTypeClass == double.class)
			{

				double[][] doubleArray = (double[][]) dataset.getData();

				for (int i = 0; i < doubleArray.length; i++)
				{
					// new vector: id, foreingKey, vector
					vectors.add(new Vector(i, 0, doubleArray[i]));
				}
				
				doubleArray = null;

			}

			else
			{
				throw new RuntimeException("Unsupported data type: " + dataTypeClass.getName());
			}
			
			hdfFile.close();

			System.out.println("[INFO] Data converted to vector objects");
			// normalize vectors
			System.out.println("[INFO] Normalizing vectors to unit size...");
			Vector.normalizeVectors(vectors);
			System.out.println("[INFO] Finished Normalization");
	
		} catch (Exception e)
		{
			System.out.print("Error occurred while reading hdf5 file: " + e);
		}

		return vectors;

	}
	
	
	public static List<Vector> fromFvecsFile(Path sourcePath) 
	{
		
        List<Vector> vectors = new ArrayList<>();
        
        try (FileChannel fc = (FileChannel.open(sourcePath, StandardOpenOption.READ))) 
        {
        	
            ByteBuffer bb = ByteBuffer.allocateDirect((int)fc.size());
            
            while (fc.read(bb) > 0);

            // Switch to read mode
            bb.flip();
            bb.order(ByteOrder.LITTLE_ENDIAN);

            while (bb.hasRemaining()) 
            {
            	
                int dimension = bb.getInt(); // Read dimension
                double[] vector = new double[dimension];
                
                for (int i = 0; i < dimension; i++) 
                {
                    vector[i] = bb.getFloat(); // Convert float to double
                }
                
                vectors.add(new Vector(vectors.size()+1, 0, vector)); // ID, ForeignKey, Vector
            }
        } catch (IOException e) 
        {
            System.err.println("Error occurred while reading .fvecs file: " + e.getMessage());
            e.printStackTrace();
        }

        // Optional: Normalize vectors
        System.out.println("[INFO] Normalizing vectors to unit size...");
        Vector.normalizeVectors(vectors);
        System.out.println("[INFO] Finished Normalization");
        System.out.println("[INFO] Vectors loaded from fvecs file successfully. Total: "+vectors.size());

        return vectors;
    }

}