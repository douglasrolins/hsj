package util;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.util.CombinatoricsUtils;

import hnswlib.Vector;

public class MetricCalculator
{

	List<Vector> vectors = new ArrayList<>();
	private List<String> results;

	public MetricCalculator(List<Vector> vectors, List<String> results)
	{
		this.vectors = vectors;
		this.results = results;
	}

	// Obtém o número total de pares positivos no dataset com base na chave estrangeira
	public int getTotalPositivesInDataset()
	{

		//# identificar quais registro possuem duplicatas e quantas
		List<Integer> duplicatas = new ArrayList<>();
		Set<Integer> unicos = new HashSet<>();
		for (Vector vector : vectors)
		{
			unicos.add(vector.foreingKey());
		}
		Map<Integer, Integer> vc = new HashMap<>();
		for (Integer id : unicos)
		{
			int count = 0;
			for (Vector vector : vectors)
			{
				if (vector.foreingKey() == id)
				{
					count++;
				}
			}
			vc.put(id, count);
		}
		for (Map.Entry<Integer, Integer> entry : vc.entrySet())
		{
			if (entry.getValue() > 1)
			{
				duplicatas.add(entry.getKey());
			}
		}

		// encontrar a quantidade total de pares positivos Combinação de(n,k)
		int total_combinacoes = 0;
		for (Integer valor : duplicatas)
		{
			int count = vc.get(valor);
			//total_combinacoes += MathUtil.comb(count, 2); // usando a classe MathUtil para calcular a combinação
			total_combinacoes += CombinatoricsUtils.binomialCoefficient(count, 2);
		}
		return total_combinacoes;
	}

	public int getTotalCombinations()
	{
		int combinations;
		combinations = (int) CombinatoricsUtils.binomialCoefficient(vectors.size(), 2);
		return combinations;
	}

	// Obtém total de positivos no resultado
	public int getTotalPositivesInResult()
	{
		int pos = 0;
		for (String r : results)
		{
			String[] elements = r.split(",");
			String third = elements[2].trim();
			String fourth = elements[3].trim();
			if (third.equals(fourth))
			{
				pos++;
			}
		}
		return pos;
	}

	public double[] getMetrics()
	{

		int tp = getTotalPositivesInResult();
		int totalResult = results.size();
		int fp = totalResult - tp;

		int totalPositivesInDataset = getTotalPositivesInDataset();
		int fn = totalPositivesInDataset - tp;

		double precisao = (double) tp / (tp + fp);
		double recall = (double) tp / (tp + fn);
		double f1 = 2 * precisao * recall / (precisao + recall);

		double[] metrics = new double[7];
		metrics[0] = precisao;
		metrics[1] = recall;
		metrics[2] = f1;
		metrics[3] = totalPositivesInDataset;
		metrics[4] = tp;
		metrics[5] = fp;
		metrics[6] = fn;

		return metrics;
	}

}