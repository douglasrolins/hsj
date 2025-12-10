import java.io.IOException;
import java.math.RoundingMode;
import java.nio.file.InvalidPathException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicLong;

import hnswlib.DistanceFunctions;
import hnswlib.SearchResult;
import hnswlib.SearchResultHops;
import hnswlib.SearchResultTotal;
import hnswlib.Vector;
import hnswlib.bruteforce.BruteForceIndex;
import hnswlib.hnsw.HnswIndex;
import util.MetricCalculator;
import util.SimException;
import util.Simcalc;

public class SimJoin
{
	private SimJoin()
	{
	}

	public static void main(String[] args) throws SimException, IOException
	{

		if (args.length < 1)
		{
			System.out.println("Usage:");
			System.out.println();
			System.out.println("  SimJoin plan=\"<algorithm>\" fileDir=\"<data directory>\" threshold=<value>");
			System.out.println("          runTimes=<int> printResult=<true|false>");
			System.out.println("          M=<int> efS=<int> efC=<int>");
			System.out.println();
			System.out.println("Parameters:");
			System.out.println("  plan        = Algorithm version to execute. Options:");
			System.out.println("                  hsj-ext  - External HSJ: operates on top of a standard HNSW index,");
			System.out.println("                             issuing multiple top-k searches until the similarity");
			System.out.println("                             threshold is reached.");
			System.out.println("                  hsj-ths  - Threshold HSJ: modified HNSW search procedure that");
			System.out.println("                             supports direct threshold-search instead of top-k.");
			System.out.println("                  hsj-inc  - Incremental HSJ: extends HSJ-Ths by interleaving");
			System.out
					.println("                             incremental index construction and threshold-based search.");
			System.out.println();
			System.out.println("  fileDir     = Path/URL of the dataset to load");
			System.out.println("  datasetHdf5     = When dataset is HDF5 select base [train, test, etc]");
			System.out.println("  threshold   = Similarity threshold (double)");
			System.out.println("  runTimes    = Number of repeated executions");
			System.out.println("  printResult = Whether to print resulting join pairs (true/false)");
			System.out.println("  M           = HNSW parameter M (max edges per node)");
			System.out.println("  ef         = HNSW search parameter efSearch");
			System.out.println("  efConstruction         = HNSW construction parameter efConstruction");
			System.exit(1);
		}

		String plan = null;
		String fileDir = null;
		String[] keyValue = null;
		String key = null;

		Map<String, String> paramsMap = new HashMap<String, String>();

		for (String param : args)
		{
			keyValue = param.split("=");
			if (keyValue.length != 2)
			{
				System.out.println("Wrong param format.");
				System.out.println("Correct param format: key=value");
				System.out.println(
						"  - There must be no white space caracter (i.e., space and tab) surrounding the delimiter \"=\".");
				System.out.println(
						"  - There must be no white space caracter in value (all valid keys do not contain white space caracters as well).");
				System.exit(1);
			}

			key = keyValue[0];

			switch (key)
			{
			case "plan":
			{
				plan = keyValue[1];
				break;
			}
			case "fileDir":
			{
				fileDir = keyValue[1];
			}

			default:
			{
				paramsMap.put(key, keyValue[1]);
			}
			}
		}

		Path sourcePath = null;
		try
		{
			sourcePath = Paths.get(fileDir);
		} catch (InvalidPathException ex)
		{
			System.out.println("Invalid path provided for fileDir: " + fileDir);
			ex.printStackTrace();
			System.exit(1);
		}

		String thsStr = paramsMap.get("threshold");
		if (thsStr == null)
		{
			System.out.println("Missing parameter: threshold");
			System.out.println("Usage example: threshold=0.85");
			System.exit(1);
		}

		double ths = 0.0;

		try
		{
			ths = Double.parseDouble(thsStr);
		} catch (NumberFormatException ex)
		{
			System.out.println("Invalid threshold format: " + thsStr);
			System.out.println("Expected a numeric value, e.g., threshold=0.85");
			System.exit(1);
		}

		if (plan == null)
		{
			System.out.println("Required parameter not specified: plan=\"query execution plan\"");
			System.exit(1);
		}

		boolean printResult = false;
		String prStr = paramsMap.get("printResult");
		if (prStr != null)
		{
			printResult = Boolean.parseBoolean(prStr);
		}

		int runTimes = 0;
		String rtStr = paramsMap.get("runTimes");
		if (rtStr != null)
		{
			try
			{
				runTimes = Integer.parseInt(rtStr);

				if (runTimes < 1)
				{
					System.out.println("Warning: \"runTimes\" param value less than 1. Using default (0).");
					runTimes = 0;
				}
			} catch (NumberFormatException e)
			{
				System.out.println("Wrong number format for \"runTimes\" param value. Using default (0).");
			}
		}

		if (!printResult && runTimes == 0)
		{
			System.out.println(
					"Neither result printing nor performance measurement were selected. Try making \"printResult=true\", or runTimes=i (i > 0), or both");
			System.exit(1);
		}

		// set dataset in hdf5 file
		String datasetHdf5 = "train"; // default
		String paramDatasetHdf5 = paramsMap.get("datasetHdf5");
		if (paramDatasetHdf5 != null)
		{
			datasetHdf5 = paramDatasetHdf5;
		}

		// HNSW parameters

		int M = 32; //default
		String mStr = paramsMap.get("M");
		if (mStr != null)
		{
			try
			{
				M = Integer.parseInt(mStr);
			} catch (Exception e)
			{
				System.out.println("M must be integer");
				System.exit(1);
			}
		}

		int efConstruction = 16; //default
		String efConstructionStr = paramsMap.get("efConstruction");
		if (efConstructionStr != null)
		{
			try
			{
				efConstruction = Integer.parseInt(efConstructionStr);
			} catch (Exception e)
			{
				System.out.println("efConstruction must be integer");
				System.exit(1);
			}
		}

		int ef = 16; //default
		String efStr = paramsMap.get("ef");
		if (efStr != null)
		{
			try
			{
				ef = Integer.parseInt(efStr);
			} catch (Exception e)
			{
				System.out.println("ef must be integer");
				System.exit(1);
			}
		}

		switch (plan)
		{

		// Baseline version, which performs multiple queries on search(from the highest layer) until the last pair is below the threshold
		// Heuristic for initial k and its increment
		// Version with parallelism (hnsw0-par -> hsj-ext)
		case "hsj-ext":
		{

			///////////
			/// read file directly and convert to vector list

			long startLoadFile = System.currentTimeMillis();

			List<Vector> vectors = new ArrayList<>();

			String typeFileDataset = sourcePath.getFileName().toString().split("\\.")[1];

			switch (typeFileDataset)
			{

			case "csv":
			{
				vectors = Vector.fromCsvFile(sourcePath);
				break;
			}

			case "hdf5":
			{
				vectors = Vector.fromHdf5File(sourcePath, datasetHdf5);
				break;
			}

			case "fvecs":
			{
				vectors = Vector.fromFvecsFile(sourcePath);
				break;
			}

			default:
				System.out.println("Unsupported file.");

			}

			long endLoadFile = System.currentTimeMillis();

			double timeLoadFile = ((endLoadFile - startLoadFile) / 1000d);

			long startIndex = System.currentTimeMillis();

			int dimension = vectors.get(0).dimensions();

			HnswIndex<Integer, double[], Vector, Double> hnswIndex = HnswIndex
					.newBuilder(dimension, DistanceFunctions.DOUBLE_COSINE_DISTANCE, vectors.size()).withM(M).withEf(ef)
					.withEfConstruction(efConstruction).build();

			try
			{
				hnswIndex.addAll(vectors);
			} catch (InterruptedException e)
			{
				e.printStackTrace();
			}

			long endIndex = System.currentTimeMillis();

			double timeIndex = ((endIndex - startIndex) / 1000d);

			//Get index memory size
			//long sizeIndex = GraphLayout.parseInstance(hnswIndex).totalSize();

			//adjust to distance
			DecimalFormat df = new DecimalFormat("0.00");
			df.setRoundingMode(RoundingMode.UP);
			//double tt = 1 - ths;
			//double threshold = Double.parseDouble(df.format(tt).replaceAll(",", "."));
			//System.out.println(ths+" - "+tt+" - "+threshold);
			//System.out.println(hnswIndex.size()+" - "+ths);

			long start;
			long end;
			double time;

			//long totalPairs = 0;

			long visitedVertices = 0;
			;
			long hops = 0;
			;
			SearchResultHops resultsHops;

			//runPerfBench
			if (runTimes > 0)
			{

				//List<SearchResult<Vector, Double>> annResults;

				//List<String> resultsFull = new ArrayList<>();
				List<String> results = new ArrayList<>();

				double[] res = new double[runTimes];

				for (int i = 0; i < runTimes; i++)
				{

					start = System.currentTimeMillis();

					// initial k - heuristic 1: k adjustable according to the ef parameter and the threshold. The lower the threshold, the higher the k.
					//int init_k = ef;
					int init_k = (int) Math.round(ef + ((1 - ths) * ef));

					try
					{
						resultsHops = hnswIndex.searchAllIndexTopK(vectors, init_k, ths);

						visitedVertices = resultsHops.getVisitedVertices();
						hops = resultsHops.getHops();
						results = resultsHops.getResults();

						//results.clear();

						/*for (String result : resultsFull) {
						    String[] parts = result.split(",");
						    int id1 = Integer.parseInt(parts[0]);
						    int id2 = Integer.parseInt(parts[1]);
						
						    if (id2 > id1) {
						        results.add(result);
						    }
						}*/

					} catch (InterruptedException e)
					{
						e.printStackTrace();
					}

					end = System.currentTimeMillis();

					res[i] = ((end - start) / 1000d);
					System.out.println(String.format("Wall-clock time of run no. %d: %.3f seconds", i, res[i]));
				}

				if (typeFileDataset.equals("hdf5") || typeFileDataset.equals("fvecs"))
				{
					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("----- Run's Information ----");
					System.out.println("Thresholds: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					time = Simcalc.calcMean(res, true);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Average wall-clock time for join: %.3f seconds.", time));
					//System.out.println("Total Pairs in Result: " + results.size());
					System.out.println("Total Pairs in Result: " + results.size());

					System.out.println("Visited Vertices: " + visitedVertices);
					System.out.println("Hops: " + hops);
				}

				else
				{
					// get foreignKey
					Map<Integer, Vector> vectorMap = new HashMap<>();
					for (Vector vector : vectors)
					{
						vectorMap.put(vector.id(), vector);
					}

					for (int i = 0; i < results.size(); i++)
					{
						String r = results.get(i);

						String[] rSplit = r.split(",");
						int idLeft = Integer.parseInt(rSplit[0]);
						int idRight = Integer.parseInt(rSplit[1]);

						// looks up the corresponding vector object in the HashMap
						Vector vectorLeft = vectorMap.get(idLeft);
						Vector vectorRight = vectorMap.get(idRight);

						// get the foreignkey of the found vector objects
						int foreignKeyLeft = vectorLeft.foreingKey();
						int foreignKeyRight = vectorRight.foreingKey();

						r = rSplit[0] + ", " + rSplit[1] + ", " + foreignKeyLeft + ", " + foreignKeyRight + ", "
								+ rSplit[2];

						results.set(i, r);
					}

					//Calculate Metrics for result
					MetricCalculator mc = new MetricCalculator(vectors, results);
					double[] metrics = mc.getMetrics();

					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("Total Positive Pairs in Dataset: (TP+FN): " + (int) metrics[3]);

					System.out.println("----- Run's Information ----");
					System.out.println("Thresholds: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					time = Simcalc.calcMean(res, true);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Average wall-clock time for join: %.3f seconds.", time));
					//System.out.println("Index memory size: " + sizeIndex/1024 + " Kb");

					System.out.println("Total Pairs in Result: " + results.size());
					System.out.println("Positive Pairs (TP): " + (int) metrics[4]);
					System.out.println("Negative Pairs (FP): " + (int) metrics[5]);
					System.out.println("Pairs Out (FN): " + (int) metrics[6]);
					System.out.printf("Precision: %.2f\n", metrics[0] * 100);
					System.out.printf("Recall: %.2f\n", metrics[1] * 100);
					System.out.printf("F1: %.2f\n", metrics[2] * 100);

					System.out.println("Visited Vertices: " + visitedVertices);
					System.out.println("Hops: " + hops);
				}

				// temporary
				//annResults = null;
				hnswIndex = null;
				vectors = null;

			}

			//PrintResult
			if (printResult)
			{
				start = System.currentTimeMillis();

				List<String> results = new ArrayList<>();
				List<String> resultsFull = new ArrayList<>();

				// initial k
				//int init_k = ef;
				int init_k = (int) Math.round(ef + ((1 - ths) * ef));

				try
				{
					resultsHops = hnswIndex.searchAllIndexTopK(vectors, init_k, ths);

					visitedVertices = resultsHops.getVisitedVertices();
					hops = resultsHops.getHops();
					resultsFull = resultsHops.getResults();

					results.clear();

					for (String result : resultsFull)
					{
						String[] parts = result.split(",");
						int id1 = Integer.parseInt(parts[0]);
						int id2 = Integer.parseInt(parts[1]);

						if (id2 > id1)
						{
							results.add(result);
						}
					}

				} catch (InterruptedException e)
				{
					e.printStackTrace();
				}

				end = System.currentTimeMillis();
				time = ((end - start) / 1000d);

				// Print results
				for (String r : results)
				{
					System.out.println(r);
				}

				if (typeFileDataset.equals("hdf5") || typeFileDataset.equals("fvecs"))
				{

					// Print results
					for (String r : results)
					{
						System.out.println(r);
					}

					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("----- Run's Information ----");
					System.out.println("Threshold: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);
					System.out.printf("Total time: %.3f seconds\n", timeLoadFile + timeIndex + time);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Wall-clock time for join operation: %.3f seconds", time));
					System.out.println("Total Pairs in Result: " + results.size());

					System.out.println("Visited Vertices: " + visitedVertices);
					System.out.println("Hops: " + hops);
				}

				else
				{

					// get foreignKey
					Map<Integer, Vector> vectorMap = new HashMap<>();
					for (Vector vector : vectors)
					{
						vectorMap.put(vector.id(), vector);
					}

					for (int i = 0; i < results.size(); i++)
					{
						String r = results.get(i);

						String[] rSplit = r.split(",");
						int idLeft = Integer.parseInt(rSplit[0]);
						int idRight = Integer.parseInt(rSplit[1]);

						// looks up the corresponding vector object in the HashMap
						Vector vectorLeft = vectorMap.get(idLeft);
						Vector vectorRight = vectorMap.get(idRight);

						// get the foreignkey of the found vector objects
						int foreignKeyLeft = vectorLeft.foreingKey();
						int foreignKeyRight = vectorRight.foreingKey();

						r = rSplit[0] + ", " + rSplit[1] + ", " + foreignKeyLeft + ", " + foreignKeyRight + ", "
								+ rSplit[2];

						results.set(i, r);

						System.out.println(r);
					}

					//Calculate Metrics for result
					MetricCalculator mc = new MetricCalculator(vectors, results);
					double[] metrics = mc.getMetrics();

					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("Total Positive Pairs in Dataset: (TP+FN): " + (int) metrics[3]);

					System.out.println("----- Run's Information ----");
					System.out.println("Threshold: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					System.out.printf("Total time: %.3f seconds\n", timeLoadFile + timeIndex + time);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Wall-clock time for join operation: %.3f seconds", time));
					//System.out.println("Index memory size: " + sizeIndex/1024 + " Kb");

					System.out.println("Total Pairs in Result: " + results.size());
					System.out.println("Positive Pairs (TP): " + (int) metrics[4]);
					System.out.println("Negative Pairs (FP): " + (int) metrics[5]);
					System.out.println("Pairs Out (FN): " + (int) metrics[6]);
					System.out.printf("Precision: %.2f\n", metrics[0] * 100);
					System.out.printf("Recall: %.2f\n", metrics[1] * 100);
					System.out.printf("F1: %.2f\n", metrics[2] * 100);

					System.out.println("Visited Vertices:  " + visitedVertices);
					System.out.println("Hops: " + hops);
				}

			}

			break;

		}

		// Baseline version, which performs multiple queries on search(from the highest layer) until the last pair is below the threshold
		// Heuristic for initial k and its increment
		// Version with parallelism and get only total pair (hnsw0-par-total -> hsj-ext-total)
		case "hsj-ext-total":
		{

			///////////
			/// read file directly and convert to vector list

			long startLoadFile = System.currentTimeMillis();

			List<Vector> vectors = new ArrayList<>();

			String typeFileDataset = sourcePath.getFileName().toString().split("\\.")[1];

			switch (typeFileDataset)
			{

			case "csv":
			{
				vectors = Vector.fromCsvFile(sourcePath);
				break;
			}

			case "hdf5":
			{
				vectors = Vector.fromHdf5File(sourcePath, datasetHdf5);
				break;
			}

			case "fvecs":
			{
				vectors = Vector.fromFvecsFile(sourcePath);
				break;
			}

			default:
				System.out.println("Unsupported file.");

			}

			long endLoadFile = System.currentTimeMillis();

			double timeLoadFile = ((endLoadFile - startLoadFile) / 1000d);

			long startIndex = System.currentTimeMillis();

			int dimension = vectors.get(0).dimensions();

			HnswIndex<Integer, double[], Vector, Double> hnswIndex = HnswIndex
					.newBuilder(dimension, DistanceFunctions.DOUBLE_COSINE_DISTANCE, vectors.size()).withM(M).withEf(ef)
					.withEfConstruction(efConstruction).build();

			try
			{
				hnswIndex.addAll(vectors);
			} catch (InterruptedException e)
			{
				e.printStackTrace();
			}

			long endIndex = System.currentTimeMillis();

			double timeIndex = ((endIndex - startIndex) / 1000d);

			//Get index memory size
			//long sizeIndex = GraphLayout.parseInstance(hnswIndex).totalSize();

			//adjust to distance
			DecimalFormat df = new DecimalFormat("0.00");
			df.setRoundingMode(RoundingMode.UP);
			//double tt = 1 - ths;
			//double threshold = Double.parseDouble(df.format(tt).replaceAll(",", "."));
			//System.out.println(ths+" - "+tt+" - "+threshold);
			//////////////////////////////////////////////////

			/////////////
			////////////

			//System.out.println(hnswIndex.size()+" - "+ths);

			long start;
			long end;
			double time;

			// hop is the number of points visited in the search route 
			long visitedVertices = 0;
			long hops = 0;
			long totalPairs = 0;
			SearchResultTotal resultsHopsTotal;

			//runPerfBench
			if (runTimes > 0)
			{

				//List<SearchResult<Vector, Double>> annResults;

				//List<String> results = new ArrayList<>();

				//int totalPairs = 0;

				double[] res = new double[runTimes];

				for (int i = 0; i < runTimes; i++)
				{

					start = System.currentTimeMillis();

					// initial k - heuristic 1: k adjustable according to the ef parameter and the threshold. The lower the threshold, the higher the k.
					//int init_k = ef;
					int init_k = (int) Math.round(ef + ((1 - ths) * ef));

					try
					{
						resultsHopsTotal = hnswIndex.searchAllIndexTopKTotal(vectors, init_k, ths);
						totalPairs = resultsHopsTotal.getTotalPairs();
						visitedVertices = resultsHopsTotal.getVisitedVertices();
						hops = resultsHopsTotal.getHops();
					} catch (InterruptedException e)
					{
						e.printStackTrace();
					}

					end = System.currentTimeMillis();

					res[i] = ((end - start) / 1000d);
					System.out.println(String.format("Wall-clock time of run no. %d: %.3f seconds", i, res[i]));
				}

				if (typeFileDataset.equals("hdf5") || typeFileDataset.equals("fvecs"))
				{
					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("----- Run's Information ----");
					System.out.println("Thresholds: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					time = Simcalc.calcMean(res, true);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Average wall-clock time for join: %.3f seconds.", time));
					System.out.println("Total Pairs in Result: " + totalPairs);

					System.out.println("Visited Vertices: " + visitedVertices);
					System.out.println("Hops: " + hops);
				}

				else
				{
					System.out.println(
							"MSG: hnsw0-par-total supported only keyless dataset! (Ex: hdf5 from ann-benchmarks");
				}

				// temporary
				//annResults = null;
				hnswIndex = null;
				vectors = null;

			}

			//PrintResult
			if (printResult)
			{
				System.out.println("MSG: hnsw0-par-total not supported printResult!");

			}

			break;

		}

		// Baseline version, which performs multiple queries on search(from the highest layer) until the last pair is below the threshold
		// Heuristic for initial k and its increment
		// Version without parallelism (hnsw0 -> hsj-ext-wopar)
		case "hsj-ext-wopar":
		{

			///////////
			/// read file directly and convert to vector list

			long startLoadFile = System.currentTimeMillis();

			List<Vector> vectors = new ArrayList<>();

			String typeFileDataset = sourcePath.getFileName().toString().split("\\.")[1];

			switch (typeFileDataset)
			{

			case "csv":
			{
				vectors = Vector.fromCsvFile(sourcePath);
				break;
			}

			case "hdf5":
			{
				vectors = Vector.fromHdf5File(sourcePath, datasetHdf5);
				break;
			}

			case "fvecs":
			{
				vectors = Vector.fromFvecsFile(sourcePath);
				break;
			}

			default:
				System.out.println("Unsupported file.");

			}

			long endLoadFile = System.currentTimeMillis();

			double timeLoadFile = ((endLoadFile - startLoadFile) / 1000d);

			long startIndex = System.currentTimeMillis();

			int dimension = vectors.get(0).dimensions();

			HnswIndex<Integer, double[], Vector, Double> hnswIndex = HnswIndex
					.newBuilder(dimension, DistanceFunctions.DOUBLE_COSINE_DISTANCE, vectors.size()).withM(M).withEf(ef)
					.withEfConstruction(efConstruction).build();

			try
			{
				hnswIndex.addAll(vectors);
			} catch (InterruptedException e)
			{
				e.printStackTrace();
			}

			long endIndex = System.currentTimeMillis();

			double timeIndex = ((endIndex - startIndex) / 1000d);

			//Get index memory size
			//long sizeIndex = GraphLayout.parseInstance(hnswIndex).totalSize();

			//adjust to distance
			DecimalFormat df = new DecimalFormat("0.00");
			df.setRoundingMode(RoundingMode.UP);
			double tt = 1 - ths;
			double threshold = Double.parseDouble(df.format(tt).replaceAll(",", "."));
			//System.out.println(ths+" - "+tt+" - "+threshold);
			//System.out.println(hnswIndex.size()+" - "+ths);

			long start;
			long end;
			double time;

			long totalPairs = 0;

			//AtomicLong revo = new AtomicLong(0);
			//AtomicLong totalRevo = new AtomicLong(0);

			AtomicLong visitedVertices = new AtomicLong(0);
			AtomicLong hops = new AtomicLong(0);

			//runPerfBench
			if (runTimes > 0)
			{

				List<SearchResult<Vector, Double>> annResults;

				List<String> results = new ArrayList<>();

				double[] res = new double[runTimes];

				for (int i = 0; i < runTimes; i++)
				{

					totalPairs = 0;

					start = System.currentTimeMillis();

					// initial k - heuristic 1: k adjustable according to the ef parameter and the threshold. The lower the threshold, the higher the k.
					//int init_k = ef;
					int init_k = (int) Math.round(ef + ((1 - ths) * ef));

					for (Vector v : vectors)
					{
						//System.out.println("------------ Vector: " + v.id() + " ---------------");

						double[] q = v.vector();

						int k = init_k;

						annResults = hnswIndex.findNearest(q, k, visitedVertices, hops);

						//revo.set(1);
						//System.out.print("- Search: " + revo.get());
						//System.out.println(", Hops: " + hops.get());

						//
						// While the most distant element is still below the threshold (considering distance) repeat the search incrementing the k
						//

						// get the similarity of the k element, which is the least similar 
						double simK = 1 - annResults.get(annResults.size() - 1).distance();

						while (simK > ths)
						{

							annResults = null;

							// Heuristic 1: Proportional increment between the similarity of the most distant element and the threshold
							k = k + (int) Math.round(((1 - ths) * k) / (1 - simK));

							// Heuristic 2: Logarithmic increment, the greater the distance from the last element to the threshold, the greater the increment.
							//k = k + (int) Math.round(k * Math.log(1 + simK / ths ) );

							annResults = hnswIndex.findNearest(q, k, visitedVertices, hops);

							//revo.incrementAndGet();

							//System.out.print("- Search: " + revo.get());
							//System.out.println(", Hops: " + hops.get());

							simK = 1 - annResults.get(annResults.size() - 1).distance();

						}

						//totalRevo.set(totalRevo.get() + revo.get());
						//System.out.println("TotalRevo: " + totalRevo.get());

						/*for (SearchResult<Vector, Double> result : annResults)
						{
							System.out.println(v.id()+" - " + result.item().id() + " - "  + result.distance() );
						}*/

						// done this way to only obtain the total of similar pairs and not the ids of each pair
						if (typeFileDataset.equals("hdf5") || typeFileDataset.equals("fvecs"))
						{
							// remove pairs that do not meet the threshold from the results
							annResults
									.removeIf(result -> result.distance() > threshold || result.item().id() <= v.id());
							totalPairs = totalPairs + annResults.size();
						} else
						{
							// temporary commented
							if (i + 1 == runTimes)
							{
								// remove pairs that do not meet the threshold from the results
								annResults.removeIf(
										result -> result.distance() > threshold || result.item().id() <= v.id());

								for (SearchResult<Vector, Double> result : annResults)
								{
									//int resultID = result.item().id();
									double distance = 1 - result.distance();
									//if (resultID > v.id() && distance > ths ) //not include repeated pairs and below the threshold
									//{
									String r = (v.id() + ", " + result.item().id() + ", "
											+ vectors.get(v.id() - 1).foreingKey() + ", " + result.item().foreingKey()
											+ ", " + distance);
									results.add(r);
									//}
								}
							}
						}

					}

					end = System.currentTimeMillis();

					res[i] = ((end - start) / 1000d);
					System.out.println(String.format("Wall-clock time of run no. %d: %.3f seconds", i, res[i]));
				}

				if (typeFileDataset.equals("hdf5") || typeFileDataset.equals("fvecs"))
				{
					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("----- Run's Information ----");
					System.out.println("Thresholds: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					time = Simcalc.calcMean(res, true);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Average wall-clock time for join: %.3f seconds.", time));
					//System.out.println("Total Pairs in Result: " + results.size());
					System.out.println("Total Pairs in Result: " + totalPairs);

					System.out.println("Visited Vertices: " + visitedVertices);
					System.out.println("Hops: " + hops);
				}

				else
				{
					//Calculate Metrics for result
					MetricCalculator mc = new MetricCalculator(vectors, results);
					double[] metrics = mc.getMetrics();

					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("Total Positive Pairs in Dataset: (TP+FN): " + (int) metrics[3]);

					System.out.println("----- Run's Information ----");
					System.out.println("Thresholds: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					time = Simcalc.calcMean(res, true);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Average wall-clock time for join: %.3f seconds.", time));
					//System.out.println("Index memory size: " + sizeIndex/1024 + " Kb");

					System.out.println("Total Pairs in Result: " + results.size());
					System.out.println("Positive Pairs (TP): " + (int) metrics[4]);
					System.out.println("Negative Pairs (FP): " + (int) metrics[5]);
					System.out.println("Pairs Out (FN): " + (int) metrics[6]);
					System.out.printf("Precision: %.2f\n", metrics[0] * 100);
					System.out.printf("Recall: %.2f\n", metrics[1] * 100);
					System.out.printf("F1: %.2f\n", metrics[2] * 100);

					System.out.println("Visited Vertices: " + visitedVertices);
					System.out.println("Hops: " + hops);
				}

				// temporary
				annResults = null;
				hnswIndex = null;
				vectors = null;

			}

			//PrintResult
			if (printResult)
			{
				start = System.currentTimeMillis();

				List<String> results = new ArrayList<>();

				// initial k
				//int init_k = ef;
				int init_k = (int) Math.round(ef + ((1 - ths) * ef));

				for (Vector v : vectors)
				{
					double[] q = v.vector();

					int k = init_k;

					List<SearchResult<Vector, Double>> annResults = hnswIndex.findNearest(q, k, visitedVertices, hops);

					//
					// While the most distant element is still below the threshold (considering distance) repeat the search incrementing the k
					//

					// get the similarity of the k element, which is the least similar 
					double simK = 1 - annResults.get(annResults.size() - 1).distance();

					while (simK > ths)
					{

						// Heuristic 1: Proportional increment between the similarity of the most distant element and the threshold
						k = k + (int) Math.round(((1 - ths) * k) / (1 - simK));

						// Heuristic 2: Logarithmic increment, the greater the distance from the last element to the threshold, the greater the increment.
						//k = k + (int) Math.round(k * Math.log(1 + simK / ths ) );

						annResults = hnswIndex.findNearest(q, k, visitedVertices, hops);

						simK = 1 - annResults.get(annResults.size() - 1).distance();

					}

					for (SearchResult<Vector, Double> result : annResults)
					{
						int resultID = result.item().id();
						double distance = 1 - result.distance();
						if (resultID > v.id() && distance > ths) //not include repeated pairs and below the threshold
						{
							String r = (v.id() + ", " + result.item().id() + ", " + vectors.get(v.id() - 1).foreingKey()
									+ ", " + result.item().foreingKey() + ", " + distance);
							results.add(r);
						}
					}
				}

				end = System.currentTimeMillis();
				time = ((end - start) / 1000d);

				//Calculate Metrics for result
				MetricCalculator mc = new MetricCalculator(vectors, results);
				double[] metrics = mc.getMetrics();

				// Print results
				for (String r : results)
				{
					System.out.println(r);
				}

				if (typeFileDataset.equals("hdf5") || typeFileDataset.equals("fvecs"))
				{
					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("----- Run's Information ----");
					System.out.println("Thresholds: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Average wall-clock time for join: %.3f seconds.", time));
					System.out.println("Total Pairs in Result: " + results.size());

					System.out.println("Visited Vertices: " + visitedVertices);
					System.out.println("Hops: " + hops);
				}

				else
				{
					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("Total Positive Pairs in Dataset: (TP+FN): " + (int) metrics[3]);

					System.out.println("----- Run's Information ----");
					System.out.println("Threshold: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					System.out.printf("Total time: %.3f seconds\n", timeLoadFile + timeIndex + time);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Wall-clock time for join operation: %.3f seconds", time));
					//System.out.println("Index memory size: " + sizeIndex/1024 + " Kb");

					System.out.println("Total Pairs in Result: " + results.size());
					System.out.println("Positive Pairs (TP): " + (int) metrics[4]);
					System.out.println("Negative Pairs (FP): " + (int) metrics[5]);
					System.out.println("Pairs Out (FN): " + (int) metrics[6]);
					System.out.printf("Precision: %.2f\n", metrics[0] * 100);
					System.out.printf("Recall: %.2f\n", metrics[1] * 100);
					System.out.printf("F1: %.2f\n", metrics[2] * 100);

					System.out.println("Visited Vertices: " + visitedVertices);
					System.out.println("Hops: " + hops);
				}

			}

			break;

		}

		// Baseline version, which performs multiple queries on search(from the highest layer) until the last pair is below the threshold.
		// With k increment of 1 by 1
		// Version without parallelism
		// (hnsw0-base -> hsj-ext-base)
		case "hsj-ext-base":
		{

			///////////hnsw0-base
			/// read file directly and convert to vector list

			long startLoadFile = System.currentTimeMillis();

			List<Vector> vectors = new ArrayList<>();

			String typeFileDataset = sourcePath.getFileName().toString().split("\\.")[1];

			switch (typeFileDataset)
			{

			case "csv":
			{
				vectors = Vector.fromCsvFile(sourcePath);
				break;
			}

			case "hdf5":
			{
				vectors = Vector.fromHdf5File(sourcePath, datasetHdf5);
				break;
			}

			case "fvecs":
			{
				vectors = Vector.fromFvecsFile(sourcePath);
				break;
			}

			default:
				System.out.println("Unsupported file.");

			}

			long endLoadFile = System.currentTimeMillis();

			double timeLoadFile = ((endLoadFile - startLoadFile) / 1000d);

			long startIndex = System.currentTimeMillis();

			int dimension = vectors.get(0).dimensions();

			HnswIndex<Integer, double[], Vector, Double> hnswIndex = HnswIndex
					.newBuilder(dimension, DistanceFunctions.DOUBLE_COSINE_DISTANCE, vectors.size()).withM(M).withEf(ef)
					.withEfConstruction(efConstruction).build();

			try
			{
				hnswIndex.addAll(vectors);
			} catch (InterruptedException e)
			{
				e.printStackTrace();
			}

			long endIndex = System.currentTimeMillis();

			double timeIndex = ((endIndex - startIndex) / 1000d);

			//Get index memory size
			//long sizeIndex = GraphLayout.parseInstance(hnswIndex).totalSize();

			//adjust to distance
			DecimalFormat df = new DecimalFormat("0.00");
			df.setRoundingMode(RoundingMode.UP);
			double tt = 1 - ths;
			double threshold = Double.parseDouble(df.format(tt).replaceAll(",", "."));
			//System.out.println(ths+" - "+tt+" - "+threshold);
			//System.out.println(hnswIndex.size()+" - "+ths);

			long start;
			long end;
			double time;

			long totalPairs = 0;

			AtomicLong visitedVertices = new AtomicLong(0);
			AtomicLong hops = new AtomicLong(0);

			//runPerfBench
			if (runTimes > 0)
			{

				List<SearchResult<Vector, Double>> annResults;

				List<String> results = new ArrayList<>();

				double[] res = new double[runTimes];

				for (int i = 0; i < runTimes; i++)
				{

					totalPairs = 0;

					start = System.currentTimeMillis();

					// initial k - heuristic 1: k adjustable according to the ef parameter and the threshold. The lower the threshold, the higher the k.
					int init_k = ef;
					//int init_k = (int) Math.round(ef + ((1-ths) * ef) );

					for (Vector v : vectors)
					{
						double[] q = v.vector();

						int k = init_k;

						annResults = hnswIndex.findNearest(q, k, visitedVertices, hops);

						//
						// While the most distant element is still below the threshold (considering distance) repeat the search incrementing the k
						//

						// get the similarity of the k element, which is the least similar 
						double simK = 1 - annResults.get(annResults.size() - 1).distance();

						while (simK > ths)
						{

							// Heuristic 1: Proportional increment between the similarity of the most distant element and the threshold
							//k = k + (int) Math.round( ((1 - ths) * k) / (1 - simK) );

							k = k + 1;

							// Heuristic 2: Logarithmic increment, the greater the distance from the last element to the threshold, the greater the increment.
							//k = k + (int) Math.round(k * Math.log(1 + simK / ths ) );

							annResults = hnswIndex.findNearest(q, k, visitedVertices, hops);

							simK = 1 - annResults.get(annResults.size() - 1).distance();

						}

						/*for (SearchResult<Vector, Double> result : annResults)
						{
							System.out.println(v.id()+" - " + result.item().id() + " - "  + result.distance() );
						}*/

						if (typeFileDataset.equals("hdf5") || typeFileDataset.equals("fvecs"))
						{
							// remove pairs that do not meet the threshold from the results
							annResults
									.removeIf(result -> result.distance() > threshold || result.item().id() <= v.id());
							totalPairs = totalPairs + annResults.size();
						} else
						{
							// temporary commented
							if (i + 1 == runTimes)
							{
								// remove pairs that do not meet the threshold from the results
								annResults.removeIf(
										result -> result.distance() > threshold || result.item().id() <= v.id());

								for (SearchResult<Vector, Double> result : annResults)
								{
									//int resultID = result.item().id();
									double distance = 1 - result.distance();
									//if (resultID > v.id() && distance > ths ) //not include repeated pairs and below the threshold
									//{
									String r = (v.id() + ", " + result.item().id() + ", "
											+ vectors.get(v.id() - 1).foreingKey() + ", " + result.item().foreingKey()
											+ ", " + distance);
									results.add(r);
									//}
								}
							}
						}

					}

					end = System.currentTimeMillis();

					res[i] = ((end - start) / 1000d);
					System.out.println(String.format("Wall-clock time of run no. %d: %.3f seconds", i, res[i]));
				}

				if (typeFileDataset.equals("hdf5") || typeFileDataset.equals("fvecs"))
				{
					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("----- Run's Information ----");
					System.out.println("Thresholds: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					time = Simcalc.calcMean(res, true);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Average wall-clock time for join: %.3f seconds.", time));
					//System.out.println("Total Pairs in Result: " + results.size());
					System.out.println("Total Pairs in Result: " + totalPairs);

					System.out.println("Visited Vertices: " + visitedVertices);
					System.out.println("Hops: " + hops);
				}

				else
				{
					//Calculate Metrics for result
					MetricCalculator mc = new MetricCalculator(vectors, results);
					double[] metrics = mc.getMetrics();

					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("Total Positive Pairs in Dataset: (TP+FN): " + (int) metrics[3]);

					System.out.println("----- Run's Information ----");
					System.out.println("Thresholds: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					time = Simcalc.calcMean(res, true);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Average wall-clock time for join: %.3f seconds.", time));
					//System.out.println("Index memory size: " + sizeIndex/1024 + " Kb");

					System.out.println("Total Pairs in Result: " + results.size());
					System.out.println("Positive Pairs (TP): " + (int) metrics[4]);
					System.out.println("Negative Pairs (FP): " + (int) metrics[5]);
					System.out.println("Pairs Out (FN): " + (int) metrics[6]);
					System.out.printf("Precision: %.2f\n", metrics[0] * 100);
					System.out.printf("Recall: %.2f\n", metrics[1] * 100);
					System.out.printf("F1: %.2f\n", metrics[2] * 100);

					System.out.println("Visited Vertices: " + visitedVertices);
					System.out.println("Hops: " + hops);
				}

				// temporary
				annResults = null;
				hnswIndex = null;
				vectors = null;

			}

			//PrintResult
			if (printResult)
			{
				start = System.currentTimeMillis();

				List<String> results = new ArrayList<>();

				// initial k
				int init_k = ef;
				//int init_k = (int) Math.round(ef + ((1-ths) * ef) );

				for (Vector v : vectors)
				{
					double[] q = v.vector();

					int k = init_k;

					List<SearchResult<Vector, Double>> annResults = hnswIndex.findNearest(q, k, visitedVertices, hops);

					//
					// While the most distant element is still below the threshold (considering distance) repeat the search incrementing the k
					//

					// get the similarity of the k element, which is the least similar 
					double simK = 1 - annResults.get(annResults.size() - 1).distance();

					while (simK > ths)
					{

						// Heuristic 1: Proportional increment between the similarity of the most distant element and the threshold
						//k = k + (int) Math.round( ((1 - ths) * k) / (1 - simK) );
						k = k + 1;

						// Heuristic 2: Logarithmic increment, the greater the distance from the last element to the threshold, the greater the increment.
						//k = k + (int) Math.round(k * Math.log(1 + simK / ths ) );

						annResults = hnswIndex.findNearest(q, k, visitedVertices, hops);

						simK = 1 - annResults.get(annResults.size() - 1).distance();

					}

					for (SearchResult<Vector, Double> result : annResults)
					{
						int resultID = result.item().id();
						double distance = 1 - result.distance();
						if (resultID > v.id() && distance > ths) //not include repeated pairs and below the threshold
						{
							String r = (v.id() + ", " + result.item().id() + ", " + vectors.get(v.id() - 1).foreingKey()
									+ ", " + result.item().foreingKey() + ", " + distance);
							results.add(r);
						}
					}
				}

				end = System.currentTimeMillis();
				time = ((end - start) / 1000d);

				//Calculate Metrics for result
				MetricCalculator mc = new MetricCalculator(vectors, results);
				double[] metrics = mc.getMetrics();

				// Print results
				for (String r : results)
				{
					System.out.println(r);
				}

				if (typeFileDataset.equals("hdf5") || typeFileDataset.equals("fvecs"))
				{
					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("----- Run's Information ----");
					System.out.println("Thresholds: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Average wall-clock time for join: %.3f seconds.", time));
					System.out.println("Total Pairs in Result: " + results.size());

					System.out.println("Visited Vertices: " + visitedVertices);
					System.out.println("Hops: " + hops);
				}

				else
				{
					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("Total Positive Pairs in Dataset: (TP+FN): " + (int) metrics[3]);

					System.out.println("----- Run's Information ----");
					System.out.println("Threshold: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					System.out.printf("Total time: %.3f seconds\n", timeLoadFile + timeIndex + time);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Wall-clock time for join operation: %.3f seconds", time));
					//System.out.println("Index memory size: " + sizeIndex/1024 + " Kb");

					System.out.println("Total Pairs in Result: " + results.size());
					System.out.println("Positive Pairs (TP): " + (int) metrics[4]);
					System.out.println("Negative Pairs (FP): " + (int) metrics[5]);
					System.out.println("Pairs Out (FN): " + (int) metrics[6]);
					System.out.printf("Precision: %.2f\n", metrics[0] * 100);
					System.out.printf("Recall: %.2f\n", metrics[1] * 100);
					System.out.printf("F1: %.2f\n", metrics[2] * 100);

					System.out.println("Visited Vertices: " + visitedVertices);
					System.out.println("Hops: " + hops);
				}

			}

			break;

		}

		// Deprecated
		// Basic version, which performs multiple queries on rangeSearch until the last pair is below the threshold
		case "hnsw1":
		{

			///////////
			/// read file directly and convert to vector list

			long startLoadFile = System.currentTimeMillis();

			List<Vector> vectors = new ArrayList<>();

			String typeFileDataset = sourcePath.getFileName().toString().split("\\.")[1];

			switch (typeFileDataset)
			{

			case "csv":
			{
				vectors = Vector.fromCsvFile(sourcePath);
				break;
			}

			case "hdf5":
			{
				vectors = Vector.fromHdf5File(sourcePath, datasetHdf5);
				break;
			}

			case "fvecs":
			{
				vectors = Vector.fromFvecsFile(sourcePath);
				break;
			}

			default:
				System.out.println("Unsupported file.");

			}

			long endLoadFile = System.currentTimeMillis();

			double timeLoadFile = ((endLoadFile - startLoadFile) / 1000d);

			long startIndex = System.currentTimeMillis();

			int dimension = vectors.get(0).dimensions();

			HnswIndex<Integer, double[], Vector, Double> hnswIndex = HnswIndex
					.newBuilder(dimension, DistanceFunctions.DOUBLE_COSINE_DISTANCE, vectors.size()).withM(M).withEf(ef)
					.withEfConstruction(efConstruction).build();

			try
			{
				hnswIndex.addAll(vectors);
			} catch (InterruptedException e)
			{
				e.printStackTrace();
			}

			long endIndex = System.currentTimeMillis();

			double timeIndex = ((endIndex - startIndex) / 1000d);

			//Get index memory size
			//long sizeIndex = GraphLayout.parseInstance(hnswIndex).totalSize();

			//adjust to distance
			DecimalFormat df = new DecimalFormat("0.00");
			df.setRoundingMode(RoundingMode.UP);
			double tt = 1 - ths;
			double threshold = Double.parseDouble(df.format(tt).replaceAll(",", "."));
			//System.out.println(ths+" - "+tt+" - "+threshold);
			//////////////////////////////////////////////////

			/////////////
			////////////

			//System.out.println(hnswIndex.size()+" - "+ths);

			long start;
			long end;
			double time;

			long totalPairs = 0;

			//runPerfBench
			if (runTimes > 0)
			{

				List<SearchResult<Vector, Double>> annResults;

				List<String> results = new ArrayList<>();

				double[] res = new double[runTimes];

				for (int i = 0; i < runTimes; i++)
				{
					totalPairs = 0;

					start = System.currentTimeMillis();

					for (Vector v : vectors)
					{
						double[] q = v.vector();
						// use rangeSearch or rangeSearch2
						annResults = hnswIndex.rangeSearch(q, threshold);

						if (typeFileDataset.equals("hdf5") || typeFileDataset.equals("fvecs"))
						{
							// remove repeated pairs
							annResults.removeIf(result -> result.item().id() <= v.id());
							totalPairs = totalPairs + annResults.size();
						} else
						{
							// temporary commented
							if (i + 1 == runTimes)
							{
								for (SearchResult<Vector, Double> result : annResults)
								{

									int resultID = result.item().id();

									if (resultID > v.id()) //not include repeated pairs
									{
										double distance = 1 - result.distance();
										//String r = (v.id() + ", " + result.item().id() + " (" + vectors.get(v.id() - 1).foreingKey() + ", "
										//		+ result.item().foreingKey() + "): [" + distance + "]");
										String r = (v.id() + ", " + result.item().id() + ", "
												+ vectors.get(v.id() - 1).foreingKey() + ", "
												+ result.item().foreingKey() + ", " + distance);
										results.add(r);
									}
								}
							}
						}

					}

					end = System.currentTimeMillis();

					res[i] = ((end - start) / 1000d);
					System.out.println(String.format("Wall-clock time of run no. %d: %.3f seconds", i, res[i]));
				}

				if (typeFileDataset.equals("hdf5") || typeFileDataset.equals("fvecs"))
				{
					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("----- Run's Information ----");
					System.out.println("Thresholds: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					time = Simcalc.calcMean(res, true);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Average wall-clock time for join: %.3f seconds.", time));
					System.out.println("Total Pairs in Result: " + totalPairs);

				}

				else
				{
					//Calculate Metrics for result
					MetricCalculator mc = new MetricCalculator(vectors, results);
					double[] metrics = mc.getMetrics();

					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("Total Positive Pairs in Dataset: (TP+FN): " + (int) metrics[3]);

					System.out.println("----- Run's Information ----");
					System.out.println("Thresholds: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					time = Simcalc.calcMean(res, true);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Average wall-clock time for join: %.3f seconds.", time));
					//System.out.println("Index memory size: " + sizeIndex/1024 + " Kb");

					System.out.println("Total Pairs in Result: " + results.size());
					System.out.println("Positive Pairs (TP): " + (int) metrics[4]);
					System.out.println("Negative Pairs (FP): " + (int) metrics[5]);
					System.out.println("Pairs Out (FN): " + (int) metrics[6]);
					System.out.printf("Precision: %.2f\n", metrics[0] * 100);
					System.out.printf("Recall: %.2f\n", metrics[1] * 100);
					System.out.printf("F1: %.2f\n", metrics[2] * 100);

				}

				// temporary
				annResults = null;
				hnswIndex = null;
				vectors = null;

			}

			//PrintResult
			if (printResult)
			{
				start = System.currentTimeMillis();

				List<String> results = new ArrayList<>();

				for (Vector v : vectors)
				{
					double[] q = v.vector();

					List<SearchResult<Vector, Double>> annResults = hnswIndex.rangeSearch(q, threshold);

					for (SearchResult<Vector, Double> result : annResults)
					{

						int resultID = result.item().id();

						if (resultID > v.id()) //not include repeated pairs
						{
							double distance = 1 - result.distance();
							//String r = (v.id() + ", " + result.item().id() + " (" + vectors.get(v.id() - 1).foreingKey() + ", "
							//		+ result.item().foreingKey() + "): [" + distance + "]");
							String r = (v.id() + ", " + result.item().id() + ", " + vectors.get(v.id() - 1).foreingKey()
									+ ", " + result.item().foreingKey() + ", " + distance);
							results.add(r);
						}
					}
				}

				end = System.currentTimeMillis();
				time = ((end - start) / 1000d);

				//Calculate Metrics for result
				MetricCalculator mc = new MetricCalculator(vectors, results);
				double[] metrics = mc.getMetrics();

				// Print results
				for (String r : results)
				{
					System.out.println(r);
				}

				if (typeFileDataset.equals("hdf5") || typeFileDataset.equals("fvecs"))
				{
					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");

					System.out.println("----- Run's Information ----");
					System.out.println("Threshold: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					System.out.printf("Total time: %.3f seconds\n", timeLoadFile + timeIndex + time);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Wall-clock time for join operation: %.3f seconds", time));
					System.out.println("Total Pairs in Result: " + results.size());
				}

				else
				{
					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("Total Positive Pairs in Dataset: (TP+FN): " + (int) metrics[3]);

					System.out.println("----- Run's Information ----");
					System.out.println("Threshold: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					System.out.printf("Total time: %.3f seconds\n", timeLoadFile + timeIndex + time);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Wall-clock time for join operation: %.3f seconds", time));
					//System.out.println("Index memory size: " + sizeIndex/1024 + " Kb");

					System.out.println("Total Pairs in Result: " + results.size());
					System.out.println("Positive Pairs (TP): " + (int) metrics[4]);
					System.out.println("Negative Pairs (FP): " + (int) metrics[5]);
					System.out.println("Pairs Out (FN): " + (int) metrics[6]);
					System.out.printf("Precision: %.2f\n", metrics[0] * 100);
					System.out.printf("Recall: %.2f\n", metrics[1] * 100);
					System.out.printf("F1: %.2f\n", metrics[2] * 100);
				}

			}

			break;

		}

		// Deprecated
		// version with changes in the HNSW library, adjusting the rangeSearchBaseLayer to search for similar pairs above the threshold
		case "hnsw2":
		{

			///////////
			/// read file directly and convert to vector list

			long startLoadFile = System.currentTimeMillis();

			List<Vector> vectors = new ArrayList<>();

			String typeFileDataset = sourcePath.getFileName().toString().split("\\.")[1];

			switch (typeFileDataset)
			{

			case "csv":
			{
				vectors = Vector.fromCsvFile(sourcePath);
				break;
			}

			case "hdf5":
			{
				vectors = Vector.fromHdf5File(sourcePath, datasetHdf5);
				break;
			}

			case "fvecs":
			{
				vectors = Vector.fromFvecsFile(sourcePath);
				break;
			}

			default:
				System.out.println("Unsupported file.");

			}

			long endLoadFile = System.currentTimeMillis();

			double timeLoadFile = ((endLoadFile - startLoadFile) / 1000d);

			long startIndex = System.currentTimeMillis();

			int dimension = vectors.get(0).dimensions();

			HnswIndex<Integer, double[], Vector, Double> hnswIndex = HnswIndex
					.newBuilder(dimension, DistanceFunctions.DOUBLE_COSINE_DISTANCE, vectors.size()).withM(M).withEf(ef)
					.withEfConstruction(efConstruction).build();

			try
			{
				hnswIndex.addAll(vectors);
			} catch (InterruptedException e)
			{
				e.printStackTrace();
			}

			long endIndex = System.currentTimeMillis();

			double timeIndex = ((endIndex - startIndex) / 1000d);

			//Get index memory size
			//long sizeIndex = GraphLayout.parseInstance(hnswIndex).totalSize();

			//adjust to distance
			DecimalFormat df = new DecimalFormat("0.00");
			df.setRoundingMode(RoundingMode.UP);
			double tt = 1 - ths;
			double threshold = Double.parseDouble(df.format(tt).replaceAll(",", "."));
			//System.out.println(ths+" - "+tt+" - "+threshold);
			//////////////////////////////////////////////////

			/////////////
			////////////

			//System.out.println(hnswIndex.size()+" - "+ths);

			long start;
			long end;
			double time;

			long totalPairs = 0;

			//runPerfBench
			if (runTimes > 0)
			{

				List<SearchResult<Vector, Double>> annResults;

				List<String> results = new ArrayList<>();

				double[] res = new double[runTimes];

				for (int i = 0; i < runTimes; i++)
				{
					totalPairs = 0;

					start = System.currentTimeMillis();

					for (Vector v : vectors)
					{
						double[] q = v.vector();
						// use rangeSearch or rangeSearch2
						annResults = hnswIndex.rangeSearch2(q, threshold);

						if (typeFileDataset.equals("hdf5") || typeFileDataset.equals("fvecs"))
						{
							// remove repeated pairs
							annResults.removeIf(result -> result.item().id() <= v.id());
							totalPairs = totalPairs + annResults.size();
						} else
						{
							if (i + 1 == runTimes)
							{
								for (SearchResult<Vector, Double> result : annResults)
								{

									int resultID = result.item().id();

									if (resultID > v.id()) //not include repeated pairs
									{
										double distance = 1 - result.distance();
										//String r = (v.id() + ", " + result.item().id() + " (" + vectors.get(v.id() - 1).foreingKey() + ", "
										//		+ result.item().foreingKey() + "): [" + distance + "]");
										String r = (v.id() + ", " + result.item().id() + ", "
												+ vectors.get(v.id() - 1).foreingKey() + ", "
												+ result.item().foreingKey() + ", " + distance);
										results.add(r);
									}
								}
							}
						}

					}

					end = System.currentTimeMillis();

					res[i] = ((end - start) / 1000d);
					System.out.println(String.format("Wall-clock time of run no. %d: %.3f seconds", i, res[i]));
				}

				if (typeFileDataset.equals("hdf5") || typeFileDataset.equals("fvecs"))
				{
					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("----- Run's Information ----");
					System.out.println("Thresholds: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					time = Simcalc.calcMean(res, true);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Average wall-clock time for join: %.3f seconds.", time));
					System.out.println("Total Pairs in Result: " + totalPairs);
				}

				else
				{
					//Calculate Metrics for result
					MetricCalculator mc = new MetricCalculator(vectors, results);
					double[] metrics = mc.getMetrics();

					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("Total Positive Pairs in Dataset: (TP+FN): " + (int) metrics[3]);

					System.out.println("----- Run's Information ----");
					System.out.println("Thresholds: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					time = Simcalc.calcMean(res, true);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Average wall-clock time for join: %.3f seconds.", time));
					//System.out.println("Index memory size: " + sizeIndex/1024 + " Kb");

					System.out.println("Total Pairs in Result: " + results.size());
					System.out.println("Positive Pairs (TP): " + (int) metrics[4]);
					System.out.println("Negative Pairs (FP): " + (int) metrics[5]);
					System.out.println("Pairs Out (FN): " + (int) metrics[6]);
					System.out.printf("Precision: %.2f\n", metrics[0] * 100);
					System.out.printf("Recall: %.2f\n", metrics[1] * 100);
					System.out.printf("F1: %.2f\n", metrics[2] * 100);
				}

				// temporary
				annResults = null;
				hnswIndex = null;
				vectors = null;

			}

			//PrintResult
			if (printResult)
			{
				start = System.currentTimeMillis();

				List<String> results = new ArrayList<>();

				for (Vector v : vectors)
				{
					double[] q = v.vector();

					List<SearchResult<Vector, Double>> annResults = hnswIndex.rangeSearch2(q, threshold);

					for (SearchResult<Vector, Double> result : annResults)
					{

						int resultID = result.item().id();

						if (resultID > v.id()) //not include repeated pairs
						{
							double distance = 1 - result.distance();
							//String r = (v.id() + ", " + result.item().id() + " (" + vectors.get(v.id() - 1).foreingKey() + ", "
							//		+ result.item().foreingKey() + "): [" + distance + "]");
							String r = (v.id() + ", " + result.item().id() + ", " + vectors.get(v.id() - 1).foreingKey()
									+ ", " + result.item().foreingKey() + ", " + distance);
							results.add(r);
						}
					}
				}

				end = System.currentTimeMillis();
				time = ((end - start) / 1000d);

				// Print results
				for (String r : results)
				{
					System.out.println(r);
				}

				if (typeFileDataset.equals("hdf5") || typeFileDataset.equals("fvecs"))
				{
					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("----- Run's Information ----");
					System.out.println("Threshold: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);
					System.out.printf("Total time: %.3f seconds\n", timeLoadFile + timeIndex + time);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Wall-clock time for join operation: %.3f seconds", time));
					System.out.println("Total Pairs in Result: " + results.size());
				}

				else
				{
					//Calculate Metrics for result
					MetricCalculator mc = new MetricCalculator(vectors, results);
					double[] metrics = mc.getMetrics();

					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("Total Positive Pairs in Dataset: (TP+FN): " + (int) metrics[3]);

					System.out.println("----- Run's Information ----");
					System.out.println("Threshold: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					System.out.printf("Total time: %.3f seconds\n", timeLoadFile + timeIndex + time);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Wall-clock time for join operation: %.3f seconds", time));
					//System.out.println("Index memory size: " + sizeIndex/1024 + " Kb");

					System.out.println("Total Pairs in Result: " + results.size());
					System.out.println("Positive Pairs (TP): " + (int) metrics[4]);
					System.out.println("Negative Pairs (FP): " + (int) metrics[5]);
					System.out.println("Pairs Out (FN): " + (int) metrics[6]);
					System.out.printf("Precision: %.2f\n", metrics[0] * 100);
					System.out.printf("Recall: %.2f\n", metrics[1] * 100);
					System.out.printf("F1: %.2f\n", metrics[2] * 100);
				}

			}

			break;

		}

		// Version that does not use search. Directly in the indexing phase of the vectors, it searches for similar pairs above the threshold
		// (hnsw3 -> hsj-inc)
		case "hsj-inc":
		{

			///////////
			/// read file directly and convert to vector list

			long startLoadFile = System.currentTimeMillis();

			List<Vector> vectors = new ArrayList<>();

			String typeFileDataset = sourcePath.getFileName().toString().split("\\.")[1];

			switch (typeFileDataset)
			{

			case "csv":
			{
				vectors = Vector.fromCsvFile(sourcePath);
				break;
			}

			case "hdf5":
			{
				vectors = Vector.fromHdf5File(sourcePath, datasetHdf5);
				break;
			}

			case "fvecs":
			{
				vectors = Vector.fromFvecsFile(sourcePath);
				break;
			}

			default:
				System.out.println("Unsupported file.");

			}

			long endLoadFile = System.currentTimeMillis();

			double timeLoadFile = ((endLoadFile - startLoadFile) / 1000d);

			//adjust to distance
			DecimalFormat df = new DecimalFormat("0.00");
			df.setRoundingMode(RoundingMode.UP);
			double tt = 1 - ths;
			double threshold = Double.parseDouble(df.format(tt).replaceAll(",", "."));

			// visitedVertices is the number of vertices visited and similarity calculations 
			// hops is the number of points visited in the search route 
			long visitedVertices = 0;
			long hops = 0;

			List<String> results = new ArrayList<>();
			SearchResultHops resultsHops;

			//runPerfBench
			if (runTimes > 0)
			{

				long start;
				long end;
				double time;

				//List<String> results = new ArrayList<>();

				int dimension = vectors.get(0).dimensions();

				double[] res = new double[runTimes];

				// Heuristic for initial k
				//int k = ef;
				int k = (int) Math.round(ef + (threshold * ef));

				for (int i = 0; i < runTimes; i++)
				{

					start = System.currentTimeMillis();

					HnswIndex<Integer, double[], Vector, Double> hnswIndex = HnswIndex
							.newBuilder(dimension, DistanceFunctions.DOUBLE_COSINE_DISTANCE, vectors.size()).withM(M)
							.withEf(k).withEfConstruction(efConstruction).build();

					try
					{
						//results = hnswIndex.addAllJoin(vectors,threshold);
						resultsHops = hnswIndex.addAllJoin(vectors, threshold);

						visitedVertices = resultsHops.getVisitedVertices();
						hops = resultsHops.getHops();
						results = resultsHops.getResults();
					} catch (InterruptedException e)
					{
						e.printStackTrace();
					}

					end = System.currentTimeMillis();

					res[i] = ((end - start) / 1000d);
					System.out.println(String.format("Wall-clock time of run no. %d: %.3f seconds", i, res[i]));

					// temporary
					hnswIndex = null;
				}

				if (typeFileDataset.equals("hdf5") || typeFileDataset.equals("fvecs"))
				{
					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("----- Run's Information ----");
					System.out.println("Thresholds: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					time = Simcalc.calcMean(res, true);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String
							.format("Average wall-clock time for join together index creation: %.3f seconds.", time));
					System.out.println(String.format("Total time: %.3f seconds.", timeLoadFile + time));
					System.out.println("Total Pairs in Result: " + results.size());

					System.out.println("Visited Vertices: " + visitedVertices);
					System.out.println("Hops: " + hops);
				}

				else

				{
					// get foreignKey
					Map<Integer, Vector> vectorMap = new HashMap<>();
					for (Vector vector : vectors)
					{
						vectorMap.put(vector.id(), vector);
					}

					for (int i = 0; i < results.size(); i++)
					{
						String r = results.get(i);

						String[] rSplit = r.split(",");
						int idLeft = Integer.parseInt(rSplit[0]);
						int idRight = Integer.parseInt(rSplit[1]);

						// looks up the corresponding vector object in the HashMap
						Vector vectorLeft = vectorMap.get(idLeft);
						Vector vectorRight = vectorMap.get(idRight);

						// get the foreignkey of the found vector objects
						int foreignKeyLeft = vectorLeft.foreingKey();
						int foreignKeyRight = vectorRight.foreingKey();

						r = rSplit[0] + ", " + rSplit[1] + ", " + foreignKeyLeft + ", " + foreignKeyRight + ", "
								+ rSplit[2];

						results.set(i, r);
					}

					//Get index memory size
					//long sizeIndex = GraphLayout.parseInstance(hnswIndex).totalSize();

					//Calculate Metrics for result
					MetricCalculator mc = new MetricCalculator(vectors, results);
					double[] metrics = mc.getMetrics();

					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("Total Positive Pairs in Dataset: (TP+FN): " + (int) metrics[3]);

					System.out.println("----- Run's Information ----");
					System.out.println("Thresholds: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					time = Simcalc.calcMean(res, true);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String
							.format("Average wall-clock time for join together index creation: %.3f seconds.", time));
					System.out.println(String.format("Total time: %.3f seconds.", timeLoadFile + time));
					//System.out.println("Index memory size: " + sizeIndex/1024 + " Kb");

					System.out.println("Total Pairs in Result: " + results.size());
					System.out.println("Positive Pairs (TP): " + (int) metrics[4]);
					System.out.println("Negative Pairs (FP): " + (int) metrics[5]);
					System.out.println("Pairs Out (FN): " + (int) metrics[6]);
					System.out.printf("Precision: %.2f\n", metrics[0] * 100);
					System.out.printf("Recall: %.2f\n", metrics[1] * 100);
					System.out.printf("F1: %.2f\n", metrics[2] * 100);

					System.out.println("Visited Vertices: " + visitedVertices);
					System.out.println("Hops: " + hops);
				}

				// temporary
				results = null;
				vectors = null;

			}

			//PrintResult
			if (printResult)
			{

				long startJoin = System.currentTimeMillis();

				int dimension = vectors.get(0).dimensions();

				HnswIndex<Integer, double[], Vector, Double> hnswIndex = HnswIndex
						.newBuilder(dimension, DistanceFunctions.DOUBLE_COSINE_DISTANCE, vectors.size()).withM(M)
						.withEf(ef).withEfConstruction(efConstruction).build();

				//List<String> results = new ArrayList<>();
				try
				{
					resultsHops = hnswIndex.addAllJoin(vectors, threshold);

					visitedVertices = resultsHops.getVisitedVertices();
					hops = resultsHops.getHops();
					results = resultsHops.getResults();
				} catch (InterruptedException e)
				{
					e.printStackTrace();
				}

				long endJoin = System.currentTimeMillis();

				double timeJoin = ((endJoin - startJoin) / 1000d);

				//Get index memory size
				//long sizeIndex = GraphLayout.parseInstance(hnswIndex).totalSize();

				if (typeFileDataset.equals("hdf5") || typeFileDataset.equals("fvecs"))
				{
					// Print results
					for (String r : results)
					{
						System.out.println(r);
					}

					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("----- Run's Information ----");
					System.out.println("Threshold: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format(
							"Wall-clock time for join operation together index creation: %.3f seconds", timeJoin));
					System.out.println(String.format("Total time: %.3f seconds.", timeLoadFile + timeJoin));
					System.out.println("Total Pairs in Result: " + results.size());

					System.out.println("Visited Vertices: " + visitedVertices);
					System.out.println("Hops: " + hops);
				}

				else
				{
					// get foreignKey
					Map<Integer, Vector> vectorMap = new HashMap<>();
					for (Vector vector : vectors)
					{
						vectorMap.put(vector.id(), vector);
					}

					for (int i = 0; i < results.size(); i++)
					{
						String r = results.get(i);

						String[] rSplit = r.split(",");
						int idLeft = Integer.parseInt(rSplit[0]);
						int idRight = Integer.parseInt(rSplit[1]);

						// looks up the corresponding vector object in the HashMap
						Vector vectorLeft = vectorMap.get(idLeft);
						Vector vectorRight = vectorMap.get(idRight);

						// get the foreignkey of the found vector objects
						int foreignKeyLeft = vectorLeft.foreingKey();
						int foreignKeyRight = vectorRight.foreingKey();

						r = rSplit[0] + ", " + rSplit[1] + ", " + foreignKeyLeft + ", " + foreignKeyRight + ", "
								+ rSplit[2];

						results.set(i, r);

						System.out.println(r);
					}

					//Calculate Metrics for result
					MetricCalculator mc = new MetricCalculator(vectors, results);
					double[] metrics = mc.getMetrics();

					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("Total Positive Pairs in Dataset: (TP+FN): " + (int) metrics[3]);

					System.out.println("----- Run's Information ----");
					System.out.println("Threshold: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					//System.out.printf("Total time: %.3f seconds\n", timeIndex + time);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format(
							"Wall-clock time for join operation together index creation: %.3f seconds", timeJoin));
					System.out.println(String.format("Total time: %.3f seconds.", timeLoadFile + timeJoin));
					//System.out.println("Index memory size: " + sizeIndex/1024 + " Kb");

					System.out.println("Total Pairs in Result: " + results.size());
					System.out.println("Positive Pairs (TP): " + (int) metrics[4]);
					System.out.println("Negative Pairs (FP): " + (int) metrics[5]);
					System.out.println("Pairs Out (FN): " + (int) metrics[6]);
					System.out.printf("Precision: %.2f\n", metrics[0] * 100);
					System.out.printf("Recall: %.2f\n", metrics[1] * 100);
					System.out.printf("F1: %.2f\n", metrics[2] * 100);

					System.out.println("Visited Vertices: " + visitedVertices);
					System.out.println("Hops: " + hops);
				}

			}

			break;

		}

		// Version that does not use search. Directly in the indexing phase of the vectors, it searches for similar pairs above the threshold
		// return only total pairs
		// (hnsw3-total -> hsj-inc-total)
		case "hsj-inc-total":
		{

			///////////
			/// read file directly and convert to vector list

			long startLoadFile = System.currentTimeMillis();

			List<Vector> vectors = new ArrayList<>();

			String typeFileDataset = sourcePath.getFileName().toString().split("\\.")[1];

			switch (typeFileDataset)
			{

			case "csv":
			{
				vectors = Vector.fromCsvFile(sourcePath);
				break;
			}

			case "hdf5":
			{
				vectors = Vector.fromHdf5File(sourcePath, datasetHdf5);
				break;
			}

			case "fvecs":
			{
				vectors = Vector.fromFvecsFile(sourcePath);
				break;
			}

			default:
				System.out.println("Unsupported file.");

			}

			long endLoadFile = System.currentTimeMillis();

			double timeLoadFile = ((endLoadFile - startLoadFile) / 1000d);

			//adjust to distance
			DecimalFormat df = new DecimalFormat("0.00");
			df.setRoundingMode(RoundingMode.UP);
			double tt = 1 - ths;
			double threshold = Double.parseDouble(df.format(tt).replaceAll(",", "."));

			// hop is the number of points visited in the search route 
			long visitedVertices = 0;
			long hops = 0;
			//long totalpairs = 0;
			SearchResultTotal resultsHopsTotal;

			//runPerfBench
			if (runTimes > 0)
			{

				long start;
				long end;
				double time;

				//List<String> results = new ArrayList<>();

				int dimension = vectors.get(0).dimensions();

				double[] res = new double[runTimes];

				// Heuristic for initial k
				//int k = ef;
				int k = (int) Math.round(ef + (threshold * ef));

				long totalPairs = 0;

				for (int i = 0; i < runTimes; i++)
				{

					start = System.currentTimeMillis();

					HnswIndex<Integer, double[], Vector, Double> hnswIndex = HnswIndex
							.newBuilder(dimension, DistanceFunctions.DOUBLE_COSINE_DISTANCE, vectors.size()).withM(M)
							.withEf(k).withEfConstruction(efConstruction).build();

					try
					{
						//totalPairs = hnswIndex.addAllJoinTotal(vectors,threshold).get();	
						resultsHopsTotal = hnswIndex.addAllJoinTotal(vectors, threshold);
						totalPairs = resultsHopsTotal.getTotalPairs();
						visitedVertices = resultsHopsTotal.getVisitedVertices();
						hops = resultsHopsTotal.getHops();

					} catch (InterruptedException e)
					{
						e.printStackTrace();
					}

					end = System.currentTimeMillis();

					res[i] = ((end - start) / 1000d);
					System.out.println(String.format("Wall-clock time of run no. %d: %.3f seconds", i, res[i]));

					// temporary
					hnswIndex = null;
				}

				if (typeFileDataset.equals("hdf5") || typeFileDataset.equals("fvecs"))
				{
					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("----- Run's Information ----");
					System.out.println("Thresholds: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					time = Simcalc.calcMean(res, true);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String
							.format("Average wall-clock time for join together index creation: %.3f seconds.", time));
					System.out.println(String.format("Total time: %.3f seconds.", timeLoadFile + time));
					System.out.println("Total Pairs in Result: " + totalPairs);

					System.out.println("Visited Vertices: " + visitedVertices);
					System.out.println("Hops: " + hops);
				}

				else

				{
					System.out
							.println("MSG: hnsw3-total supported only keyless dataset! (Ex: hdf5 from ann-benchmarks");
				}

				// temporary
				//results = null;
				vectors = null;

			}

			//PrintResult
			if (printResult)
			{
				System.out.println("MSG: hnsw3-total not supported printResult!");
			}

			break;

		}

		// version that first indexes and then searches with rangeSearch in parallel
		// (hnsw4 -> hsj-ths)
		case "hsj-ths":
		{

			///////////
			/// read file directly and convert to vector list

			long startLoadFile = System.currentTimeMillis();

			List<Vector> vectors = new ArrayList<>();

			String typeFileDataset = sourcePath.getFileName().toString().split("\\.")[1];

			switch (typeFileDataset)
			{

			case "csv":
			{
				vectors = Vector.fromCsvFile(sourcePath);
				break;
			}

			case "hdf5":
			{
				vectors = Vector.fromHdf5File(sourcePath, datasetHdf5);
				break;
			}

			case "fvecs":
			{
				vectors = Vector.fromFvecsFile(sourcePath);
				break;
			}

			default:
				System.out.println("Unsupported file.");

			}

			long endLoadFile = System.currentTimeMillis();

			double timeLoadFile = ((endLoadFile - startLoadFile) / 1000d);

			long startIndex = System.currentTimeMillis();

			int dimension = vectors.get(0).dimensions();

			HnswIndex<Integer, double[], Vector, Double> hnswIndex = HnswIndex
					.newBuilder(dimension, DistanceFunctions.DOUBLE_COSINE_DISTANCE, vectors.size()).withM(M).withEf(ef)
					.withEfConstruction(efConstruction).build();

			try
			{
				hnswIndex.addAll(vectors);
			} catch (InterruptedException e)
			{
				e.printStackTrace();
			}

			long endIndex = System.currentTimeMillis();

			double timeIndex = ((endIndex - startIndex) / 1000d);

			//Get index memory size
			//long sizeIndex = GraphLayout.parseInstance(hnswIndex).totalSize();

			//adjust to distance
			DecimalFormat df = new DecimalFormat("0.00");
			df.setRoundingMode(RoundingMode.UP);
			double tt = 1 - ths;
			double threshold = Double.parseDouble(df.format(tt).replaceAll(",", "."));
			//System.out.println(ths+" - "+tt+" - "+threshold);
			//////////////////////////////////////////////////

			/////////////
			////////////

			//System.out.println(hnswIndex.size()+" - "+ths);

			long start;
			long end;
			double time;

			// hop is the number of points visited in the search route 
			long visitedVertices = 0;
			long hops = 0;
			List<String> results = new ArrayList<>();
			SearchResultHops resultsHops;

			//runPerfBench
			if (runTimes > 0)
			{

				//List<SearchResult<Vector, Double>> annResults;

				//List<String> results = new ArrayList<>();

				double[] res = new double[runTimes];

				for (int i = 0; i < runTimes; i++)
				{

					start = System.currentTimeMillis();

					try
					{
						resultsHops = hnswIndex.searchAllIndex(vectors, threshold);

						visitedVertices = resultsHops.getVisitedVertices();
						hops = resultsHops.getHops();
						results = resultsHops.getResults();
					} catch (InterruptedException e)
					{
						e.printStackTrace();
					}

					end = System.currentTimeMillis();

					res[i] = ((end - start) / 1000d);
					System.out.println(String.format("Wall-clock time of run no. %d: %.3f seconds", i, res[i]));
				}

				if (typeFileDataset.equals("hdf5") || typeFileDataset.equals("fvecs"))
				{
					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("----- Run's Information ----");
					System.out.println("Thresholds: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					time = Simcalc.calcMean(res, true);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Average wall-clock time for join: %.3f seconds.", time));
					System.out.println("Total Pairs in Result: " + results.size());

					System.out.println("Visited Vertices: " + visitedVertices);
					System.out.println("Hops: " + hops);

				}

				else
				{

					// get foreignKey
					Map<Integer, Vector> vectorMap = new HashMap<>();
					for (Vector vector : vectors)
					{
						vectorMap.put(vector.id(), vector);
					}

					for (int i = 0; i < results.size(); i++)
					{
						String r = results.get(i);

						String[] rSplit = r.split(",");
						int idLeft = Integer.parseInt(rSplit[0]);
						int idRight = Integer.parseInt(rSplit[1]);

						// looks up the corresponding vector object in the HashMap
						Vector vectorLeft = vectorMap.get(idLeft);
						Vector vectorRight = vectorMap.get(idRight);

						// get the foreignkey of the found vector objects
						int foreignKeyLeft = vectorLeft.foreingKey();
						int foreignKeyRight = vectorRight.foreingKey();

						r = rSplit[0] + ", " + rSplit[1] + ", " + foreignKeyLeft + ", " + foreignKeyRight + ", "
								+ rSplit[2];

						results.set(i, r);
					}

					//Calculate Metrics for result
					MetricCalculator mc = new MetricCalculator(vectors, results);
					double[] metrics = mc.getMetrics();

					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("Total Positive Pairs in Dataset: (TP+FN): " + (int) metrics[3]);

					System.out.println("----- Run's Information ----");
					System.out.println("Thresholds: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					time = Simcalc.calcMean(res, true);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Average wall-clock time for join: %.3f seconds.", time));
					//System.out.println("Index memory size: " + sizeIndex/1024 + " Kb");

					System.out.println("Total Pairs in Result: " + results.size());
					System.out.println("Positive Pairs (TP): " + (int) metrics[4]);
					System.out.println("Negative Pairs (FP): " + (int) metrics[5]);
					System.out.println("Pairs Out (FN): " + (int) metrics[6]);
					System.out.printf("Precision: %.2f\n", metrics[0] * 100);
					System.out.printf("Recall: %.2f\n", metrics[1] * 100);
					System.out.printf("F1: %.2f\n", metrics[2] * 100);

					System.out.println("Visited Vertices: " + visitedVertices);
					System.out.println("Hops: " + hops);
				}

				// temporary
				//annResults = null;
				hnswIndex = null;
				vectors = null;

			}

			//PrintResult
			if (printResult)
			{
				start = System.currentTimeMillis();

				//List<String> results = new ArrayList<>();
				try
				{
					resultsHops = hnswIndex.searchAllIndex(vectors, threshold);

					visitedVertices = resultsHops.getVisitedVertices();
					results = resultsHops.getResults();
					hops = resultsHops.getHops();
				} catch (InterruptedException e)
				{
					e.printStackTrace();
				}

				end = System.currentTimeMillis();
				time = ((end - start) / 1000d);

				if (typeFileDataset.equals("hdf5") || typeFileDataset.equals("fvecs"))
				{

					// Print results
					for (String r : results)
					{
						System.out.println(r);
					}

					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("----- Run's Information ----");
					System.out.println("Threshold: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);
					System.out.printf("Total time: %.3f seconds\n", timeLoadFile + timeIndex + time);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Wall-clock time for join operation: %.3f seconds", time));
					System.out.println("Total Pairs in Result: " + results.size());

					System.out.println("Visited Vertices: " + visitedVertices);
					System.out.println("Hops: " + hops);
				}

				else
				{

					// get foreignKey
					Map<Integer, Vector> vectorMap = new HashMap<>();
					for (Vector vector : vectors)
					{
						vectorMap.put(vector.id(), vector);
					}

					for (int i = 0; i < results.size(); i++)
					{
						String r = results.get(i);

						String[] rSplit = r.split(",");
						int idLeft = Integer.parseInt(rSplit[0]);
						int idRight = Integer.parseInt(rSplit[1]);

						// looks up the corresponding vector object in the HashMap
						Vector vectorLeft = vectorMap.get(idLeft);
						Vector vectorRight = vectorMap.get(idRight);

						// get the foreignkey of the found vector objects
						int foreignKeyLeft = vectorLeft.foreingKey();
						int foreignKeyRight = vectorRight.foreingKey();

						r = rSplit[0] + ", " + rSplit[1] + ", " + foreignKeyLeft + ", " + foreignKeyRight + ", "
								+ rSplit[2];

						results.set(i, r);

						System.out.println(r);
					}

					//Calculate Metrics for result
					MetricCalculator mc = new MetricCalculator(vectors, results);
					double[] metrics = mc.getMetrics();

					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("Total Positive Pairs in Dataset: (TP+FN): " + (int) metrics[3]);

					System.out.println("----- Run's Information ----");
					System.out.println("Threshold: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					System.out.printf("Total time: %.3f seconds\n", timeLoadFile + timeIndex + time);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Wall-clock time for join operation: %.3f seconds", time));
					//System.out.println("Index memory size: " + sizeIndex/1024 + " Kb");

					System.out.println("Total Pairs in Result: " + results.size());
					System.out.println("Positive Pairs (TP): " + (int) metrics[4]);
					System.out.println("Negative Pairs (FP): " + (int) metrics[5]);
					System.out.println("Pairs Out (FN): " + (int) metrics[6]);
					System.out.printf("Precision: %.2f\n", metrics[0] * 100);
					System.out.printf("Recall: %.2f\n", metrics[1] * 100);
					System.out.printf("F1: %.2f\n", metrics[2] * 100);

					System.out.println("Visited Vertices:  " + visitedVertices);
					System.out.println("Hops: " + hops);
				}

			}

			break;

		}

		// version that first indexes and then searches with rangeSearch in parallel
		// return only total pairs
		// (hnsw4-total -> hsj-ths-total)
		case "hsj-ths-total":
		{

			///////////
			/// read file directly and convert to vector list

			long startLoadFile = System.currentTimeMillis();

			List<Vector> vectors = new ArrayList<>();

			String typeFileDataset = sourcePath.getFileName().toString().split("\\.")[1];

			switch (typeFileDataset)
			{

			case "csv":
			{
				vectors = Vector.fromCsvFile(sourcePath);
				break;
			}

			case "hdf5":
			{
				vectors = Vector.fromHdf5File(sourcePath, datasetHdf5);
				break;
			}

			case "fvecs":
			{
				vectors = Vector.fromFvecsFile(sourcePath);
				break;
			}

			default:
				System.out.println("Unsupported file.");

			}

			long endLoadFile = System.currentTimeMillis();

			double timeLoadFile = ((endLoadFile - startLoadFile) / 1000d);

			long startIndex = System.currentTimeMillis();

			int dimension = vectors.get(0).dimensions();

			HnswIndex<Integer, double[], Vector, Double> hnswIndex = HnswIndex
					.newBuilder(dimension, DistanceFunctions.DOUBLE_COSINE_DISTANCE, vectors.size()).withM(M).withEf(ef)
					.withEfConstruction(efConstruction).build();

			try
			{
				hnswIndex.addAll(vectors);
			} catch (InterruptedException e)
			{
				e.printStackTrace();
			}

			long endIndex = System.currentTimeMillis();

			double timeIndex = ((endIndex - startIndex) / 1000d);

			//Get index memory size
			//long sizeIndex = GraphLayout.parseInstance(hnswIndex).totalSize();

			//adjust to distance
			DecimalFormat df = new DecimalFormat("0.00");
			df.setRoundingMode(RoundingMode.UP);
			double tt = 1 - ths;
			double threshold = Double.parseDouble(df.format(tt).replaceAll(",", "."));
			//System.out.println(ths+" - "+tt+" - "+threshold);
			//////////////////////////////////////////////////

			/////////////
			////////////

			//System.out.println(hnswIndex.size()+" - "+ths);

			long start;
			long end;
			double time;

			// hop is the number of points visited in the search route 
			long visitedVertices = 0;
			long hops = 0;
			long totalPairs = 0;
			SearchResultTotal resultsHopsTotal;

			//runPerfBench
			if (runTimes > 0)
			{

				//List<SearchResult<Vector, Double>> annResults;

				//List<String> results = new ArrayList<>();

				//int totalPairs = 0;

				double[] res = new double[runTimes];

				for (int i = 0; i < runTimes; i++)
				{

					start = System.currentTimeMillis();

					try
					{
						resultsHopsTotal = hnswIndex.searchAllIndexTotal(vectors, threshold);
						totalPairs = resultsHopsTotal.getTotalPairs();
						visitedVertices = resultsHopsTotal.getVisitedVertices();
						hops = resultsHopsTotal.getHops();
					} catch (InterruptedException e)
					{
						e.printStackTrace();
					}

					end = System.currentTimeMillis();

					res[i] = ((end - start) / 1000d);
					System.out.println(String.format("Wall-clock time of run no. %d: %.3f seconds", i, res[i]));
				}

				if (typeFileDataset.equals("hdf5") || typeFileDataset.equals("fvecs"))
				{
					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("----- Run's Information ----");
					System.out.println("Thresholds: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					time = Simcalc.calcMean(res, true);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Average wall-clock time for join: %.3f seconds.", time));
					System.out.println("Total Pairs in Result: " + totalPairs);

					System.out.println("Visited Vertices: " + visitedVertices);
					System.out.println("Hops: " + hops);
				}

				else
				{
					System.out
							.println("MSG: hnsw4-total supported only keyless dataset! (Ex: hdf5 from ann-benchmarks");
				}

				// temporary
				//annResults = null;
				hnswIndex = null;
				vectors = null;

			}

			//PrintResult
			if (printResult)
			{
				System.out.println("MSG: hnsw4-total not supported printResult!");

			}

			break;

		}

		// brute force version with hnswlib
		case "hnsw_bf":
		{

			///////////
			/// read file directly and convert to vector list

			long startLoadFile = System.currentTimeMillis();

			List<Vector> vectors = new ArrayList<>();

			String typeFileDataset = sourcePath.getFileName().toString().split("\\.")[1];

			switch (typeFileDataset)
			{

			case "csv":
			{
				vectors = Vector.fromCsvFile(sourcePath);
				break;
			}

			case "hdf5":
			{
				vectors = Vector.fromHdf5File(sourcePath, datasetHdf5);
				break;
			}

			case "fvecs":
			{
				vectors = Vector.fromFvecsFile(sourcePath);
				break;
			}

			default:
				System.out.println("Unsupported file.");

			}

			long endLoadFile = System.currentTimeMillis();

			double timeLoadFile = ((endLoadFile - startLoadFile) / 1000d);

			long startIndex = System.currentTimeMillis();

			int dimension = vectors.get(0).dimensions();

			//HnswIndex<Integer, double[], Vector, Double> hnswIndex = HnswIndex
			//		.newBuilder(dimension, DistanceFunctions.DOUBLE_COSINE_DISTANCE, vectors.size()).withM(M).withEf(ef)
			//		.withEfConstruction(efConstruction).build();

			//Index<Integer, double[], Vector, Double> groundTruthIndex = hnswIndex.asExactIndex();

			BruteForceIndex<Integer, double[], Vector, Double> bfIndex = BruteForceIndex
					.newBuilder(dimension, DistanceFunctions.DOUBLE_COSINE_DISTANCE).build();

			try
			{
				//hnswIndex.addAll(vectors);
				//groundTruthIndex.addAll(vectors);
				bfIndex.addAll(vectors);

			} catch (InterruptedException e)
			{
				e.printStackTrace();
			}

			long endIndex = System.currentTimeMillis();

			double timeIndex = ((endIndex - startIndex) / 1000d);

			//Get index memory size
			//long sizeIndex = GraphLayout.parseInstance(hnswIndex).totalSize();

			//adjust to distance
			DecimalFormat df = new DecimalFormat("0.00");
			df.setRoundingMode(RoundingMode.UP);
			double tt = 1 - ths;
			double threshold = Double.parseDouble(df.format(tt).replaceAll(",", "."));
			//System.out.println(ths+" - "+tt+" - "+threshold);
			//////////////////////////////////////////////////

			/////////////
			////////////

			//System.out.println(hnswIndex.size()+" - "+ths);

			long start;
			long end;
			double time;

			//runPerfBench
			if (runTimes > 0)
			{

				List<SearchResult<Vector, Double>> annResults;

				List<String> results = new ArrayList<>();

				double[] res = new double[runTimes];

				long totalPairs;

				AtomicLong visitedVertices = new AtomicLong();
				AtomicLong hops = new AtomicLong();

				for (int i = 0; i < runTimes; i++)
				{

					totalPairs = 0;

					start = System.currentTimeMillis();

					// initial k - heuristic 1: k adjustable according to the ef parameter and the threshold. The lower the threshold, the higher the k.
					//int init_k = ef;
					int init_k = (int) Math.round(ef + ((1 - ths) * ef));

					for (Vector v : vectors)
					{
						double[] q = v.vector();

						int k = init_k;

						//annResults = hnswIndex.findNearest(q, k);
						annResults = bfIndex.findNearest(q, k, visitedVertices, hops);

						//
						// While the most distant element is still below the threshold (considering distance) repeat the search incrementing the k
						//

						// get the similarity of the k element, which is the least similar 
						double simK = 1 - annResults.get(annResults.size() - 1).distance();

						while (simK > ths)
						{

							// Heuristic 1: Proportional increment between the similarity of the most distant element and the threshold
							k = k + (int) Math.round(((1 - ths) * k) / (1 - simK));

							// Heuristic 2: Logarithmic increment, the greater the distance from the last element to the threshold, the greater the increment.
							//k = k + (int) Math.round(k * Math.log(1 + simK / ths ) );

							annResults = bfIndex.findNearest(q, k, visitedVertices, hops);

							simK = 1 - annResults.get(annResults.size() - 1).distance();

						}

						if (typeFileDataset.equals("hdf5") || typeFileDataset.equals("fvecs"))
						{
							// remove pairs that do not meet the threshold from the results
							annResults
									.removeIf(result -> result.distance() > threshold || result.item().id() <= v.id());
							totalPairs = totalPairs + annResults.size();
						} else
						{
							// temporary commented
							if (i + 1 == runTimes)
							{
								// remove pairs that do not meet the threshold from the results
								annResults.removeIf(
										result -> result.distance() > threshold || result.item().id() <= v.id());

								for (SearchResult<Vector, Double> result : annResults)
								{
									//int resultID = result.item().id();
									double distance = 1 - result.distance();
									//if (resultID > v.id() && distance > ths ) //not include repeated pairs and below the threshold
									//{
									String r = (v.id() + ", " + result.item().id() + ", "
											+ vectors.get(v.id() - 1).foreingKey() + ", " + result.item().foreingKey()
											+ ", " + distance);
									results.add(r);
									//}
								}
							}
						}

					}

					end = System.currentTimeMillis();

					res[i] = ((end - start) / 1000d);
					System.out.println(String.format("Wall-clock time of run no. %d: %.3f seconds", i, res[i]));
				}

				if (typeFileDataset.equals("hdf5") || typeFileDataset.equals("fvecs"))
				{
					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("----- Run's Information ----");
					System.out.println("Thresholds: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					time = Simcalc.calcMean(res, true);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Average wall-clock time for join: %.3f seconds.", time));
					System.out.println("Total Pairs in Result: " + results.size());

					System.out.println("Visited Vertices: " + visitedVertices);
					System.out.println("Hops: " + hops);
				}

				else
				{

					//Calculate Metrics for result
					MetricCalculator mc = new MetricCalculator(vectors, results);
					double[] metrics = mc.getMetrics();

					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("Total Positive Pairs in Dataset: (TP+FN): " + (int) metrics[3]);

					System.out.println("----- Run's Information ----");
					System.out.println("Thresholds: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					time = Simcalc.calcMean(res, true);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Average wall-clock time for join: %.3f seconds.", time));
					//System.out.println("Index memory size: " + sizeIndex/1024 + " Kb");

					System.out.println("Total Pairs in Result: " + results.size());
					System.out.println("Positive Pairs (TP): " + (int) metrics[4]);
					System.out.println("Negative Pairs (FP): " + (int) metrics[5]);
					System.out.println("Pairs Out (FN): " + (int) metrics[6]);
					System.out.printf("Precision: %.2f\n", metrics[0] * 100);
					System.out.printf("Recall: %.2f\n", metrics[1] * 100);
					System.out.printf("F1: %.2f\n", metrics[2] * 100);

					System.out.println("Visited Vertices: " + visitedVertices);
					System.out.println("Hops: " + hops);
				}

				// temporary
				annResults = null;
				bfIndex = null;
				vectors = null;

			}

			//PrintResult
			if (printResult)
			{
				start = System.currentTimeMillis();

				AtomicLong visitedVertices = new AtomicLong();
				AtomicLong hops = new AtomicLong();

				List<String> results = new ArrayList<>();

				// initial k
				//int init_k = ef;
				int init_k = (int) Math.round(ef + ((1 - ths) * ef));

				for (Vector v : vectors)
				{
					double[] q = v.vector();

					int k = init_k;

					//List<SearchResult<Vector, Double>> annResults = hnswIndex.findNearest(q, k);	
					List<SearchResult<Vector, Double>> annResults = bfIndex.findNearest(q, k, visitedVertices, hops);

					//
					// While the most distant element is still below the threshold (considering distance) repeat the search incrementing the k
					//

					// get the similarity of the k element, which is the least similar 
					double simK = 1 - annResults.get(annResults.size() - 1).distance();

					while (simK > ths)
					{

						// Heuristic 1: Proportional increment between the similarity of the most distant element and the threshold
						k = k + (int) Math.round(((1 - ths) * k) / (1 - simK));

						// Heuristic 2: Logarithmic increment, the greater the distance from the last element to the threshold, the greater the increment.
						//k = k + (int) Math.round(k * Math.log(1 + simK / ths ) );

						annResults = bfIndex.findNearest(q, k, visitedVertices, hops);

						simK = 1 - annResults.get(annResults.size() - 1).distance();

					}

					for (SearchResult<Vector, Double> result : annResults)
					{
						int resultID = result.item().id();
						double distance = 1 - result.distance();
						if (resultID > v.id() && distance > ths) //not include repeated pairs and below the threshold
						{
							String r = (v.id() + ", " + result.item().id() + ", " + vectors.get(v.id() - 1).foreingKey()
									+ ", " + result.item().foreingKey() + ", " + distance);
							results.add(r);
						}
					}
				}

				end = System.currentTimeMillis();
				time = ((end - start) / 1000d);

				//Calculate Metrics for result
				MetricCalculator mc = new MetricCalculator(vectors, results);
				double[] metrics = mc.getMetrics();

				// Print results
				for (String r : results)
				{
					System.out.println(r);
				}

				if (typeFileDataset.equals("hdf5") || typeFileDataset.equals("fvecs"))
				{
					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("----- Run's Information ----");
					System.out.println("Thresholds: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Average wall-clock time for join: %.3f seconds.", time));
					System.out.println("Total Pairs in Result: " + results.size());

					System.out.println("Visited Vertices: " + visitedVertices);
					System.out.println("Hops: " + hops);
				}

				else
				{
					System.out.println("----- Dataset's Information ----");
					System.out.println("File: " + sourcePath);
					System.out.println("Total vectors and dimensions: (" + vectors.size() + "," + dimension + ")");
					System.out.println("Total Positive Pairs in Dataset: (TP+FN): " + (int) metrics[3]);

					System.out.println("----- Run's Information ----");
					System.out.println("Threshold: " + ths);
					System.out.println("Parameters: Index=" + plan + ", M=" + M + ", EfConstruction=" + efConstruction
							+ ", EfSearch=" + ef);

					System.out.printf("Total time: %.3f seconds\n", timeLoadFile + timeIndex + time);
					System.out.println(String.format("Wall-clock time for load file: %.3f seconds", timeLoadFile));
					System.out.println(String.format("Wall-clock time for create index: %.3f seconds", timeIndex));
					System.out.println(String.format("Wall-clock time for join operation: %.3f seconds", time));
					//System.out.println("Index memory size: " + sizeIndex/1024 + " Kb");

					System.out.println("Total Pairs in Result: " + results.size());
					System.out.println("Positive Pairs (TP): " + (int) metrics[4]);
					System.out.println("Negative Pairs (FP): " + (int) metrics[5]);
					System.out.println("Pairs Out (FN): " + (int) metrics[6]);
					System.out.printf("Precision: %.2f\n", metrics[0] * 100);
					System.out.printf("Recall: %.2f\n", metrics[1] * 100);
					System.out.printf("F1: %.2f\n", metrics[2] * 100);

					System.out.println("Visited Vertices: " + visitedVertices);
					System.out.println("Hops: " + hops);
				}

			}

			break;

		}

		default:
		{
			System.out.println(String.format("Unknown plan: %s", plan));
		}
		}
		System.exit(0);
	}

}