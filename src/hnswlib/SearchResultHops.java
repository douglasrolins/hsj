package hnswlib;

import java.util.List;

public class SearchResultHops {
    private long visitedVertices;
    private long hops;
    private List<String> results;

    public SearchResultHops(long visitedVertices, long hops, List<String> results) {
        this.visitedVertices = visitedVertices;
        this.hops = hops;
        this.results = results;
    }

    public long getVisitedVertices() {
        return visitedVertices;
    }
    
    public long getHops() {
    	return hops;
    }

    public List<String> getResults() {
        return results;
    }
}