package hnswlib;

//Class to store the total number of pairs returned in the query and the similarity value of the least similar pair

public class SearchResultTotal 
{
	
    private long totalPairs;
    private double minSimilarity;
    private long visitedVertices;
    private long hops;

    public SearchResultTotal(long totalPairs, double minSimilarity, long visitedVertices, long hops) 
    {
        this.totalPairs = totalPairs;
        this.minSimilarity = minSimilarity;
        this.visitedVertices = visitedVertices;
        this.hops = hops;
    }

    public long getTotalPairs() 
    {
        return totalPairs;
    }

    public double getMinSimilarity() 
    {
        return minSimilarity;
    }
    
    public long getVisitedVertices() {
    	return visitedVertices;
    }
    
    public long getHops() {
    	return hops;
    }
}