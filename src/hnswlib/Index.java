package hnswlib;

import hnswlib.util.NamedThreadFactory;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Collectors;

/**
 * K-nearest neighbors search index.
 *
 * @param <TId> Type of the external identifier of an item
 * @param <TVector> Type of the vector to perform distance calculation on
 * @param <TItem> Type of items stored in the index
 * @param <TDistance> Type of distance between items (expect any numeric type: float, double, int, ..)
 *
 * @see <a href="https://en.wikipedia.org/wiki/K-nearest_neighbors_algorithm">k-nearest neighbors algorithm</a>
 */
public interface Index<TId, TVector, TItem extends Item<TId, TVector>, TDistance> extends Serializable {

    /**
     * By default after indexing this many items progress will be reported to registered progress listeners.
     */
    int DEFAULT_PROGRESS_UPDATE_INTERVAL = 100_000;

    /**
     * Add a new item to the index. If an item with the same identifier already exists in the index then :
     *
     * If deletes are disabled on this index the method will return false and the item will not be updated.
     *
     * If deletes are enabled and the version of the item has is higher version than that of the item currently stored
     * in the index the old item will be removed and the new item added, otherwise this method will return false and the
     * item will not be updated.
     *
     * @param item the item to add to the index
     *
     * @return true if the item was added to the index
     * @throws IllegalArgumentException thrown when the item has the wrong dimensionality
     */
    boolean add(TItem item);
    
    /**
     * search a item in the the index and, in parallel, search similar items. 
     * If an item with the same identifier already exists in the index then :
     *
     * If deletes are disabled on this index the method will return false and the item will not be updated.
     *
     * If deletes are enabled and the version of the item has is higher version than that of the item currently stored
     * in the index the old item will be removed and the new item added, otherwise this method will return false and the
     * item will not be updated.
     *
     * @param item the item to add to the index
     * @param threshold the threshold used in the search
     *
     * @return true if the item was added to the index
     * @throws IllegalArgumentException thrown when the item has the wrong dimensionality
     */
    List<String> searchParallel(TItem item, TDistance distance, AtomicLong visitedVertices, AtomicLong hops);
    
    /**
     * Search a single item in the index and, in parallel, retrieve the top-k most similar items.
     * If an item with the same identifier already exists in the index, then:
     *
     * If deletes are disabled on this index the method will return false and the item will not be updated.
     *
     * If deletes are enabled and the version of the item is higher than that of the item currently stored
     * in the index, the old item will be removed and the new item added. Otherwise, this method will return false and the
     * item will not be updated.
     *
     * @param item the item to add to the index
     * @param k the number of top similar items to retrieve
     * @param visitedVertices output parameter to store the number of visited vertices during the search
     * @param hops output parameter to store the number of hops performed during the search
     *
     * @return list of similar item identifiers
     * @throws IllegalArgumentException thrown when the item has the wrong dimensionality
     */
    List<String> searchParallelTopK(TItem item, int k, AtomicLong visitedVertices, AtomicLong hops);

    
    /**
     * Add a new item to the index and, in parallel, search similar items. 
     * If an item with the same identifier already exists in the index then :
     *
     * If deletes are disabled on this index the method will return false and the item will not be updated.
     *
     * If deletes are enabled and the version of the item has is higher version than that of the item currently stored
     * in the index the old item will be removed and the new item added, otherwise this method will return false and the
     * item will not be updated.
     *
     * @param item the item to add to the index
     * @param threshold the threshold used in the search
     *
     * @return true if the item was added to the index
     * @throws IllegalArgumentException thrown when the item has the wrong dimensionality
     */
    List<String> addJoin(TItem item, TDistance distance, AtomicLong visitedVertices, AtomicLong hops);

    /**
     * Removes an item from the index. If the index does not support deletes or an item with the same identifier exists
     * in the index with a higher version number, then this method will return false and the item will not be removed.
     *
     * @param id unique identifier or the item to remove
     * @param version version of the delete. If your items don't override version  use 0
     * @return {@code true} if an item was removed from the index. In case the index does not support removals this will
     *                      always be false
     */
    boolean remove(TId id, long version);

    /**
     * Check if an item is contained in this index
     *
     * @param id unique identifier of the item
     * @return true if an item is contained in this index, false otherwise
     */
    default boolean contains(TId id) {
        return get(id).isPresent();
    }

    /**
     * Add multiple items to the index
     *
     * @param items the items to add to the index
     * @throws InterruptedException thrown when the thread doing the indexing is interrupted
     */
    default void addAll(Collection<TItem> items) throws InterruptedException {
        addAll(items, NullProgressListener.INSTANCE);
    }

    /**
     * Add multiple items to the index. Reports progress to the passed in implementation of {@link ProgressListener}
     * every {@link Index#DEFAULT_PROGRESS_UPDATE_INTERVAL} elements indexed.
     *
     * @param items the items to add to the index
     * @param listener listener to report progress to
     * @throws InterruptedException thrown when the thread doing the indexing is interrupted
     */
    default void addAll(Collection<TItem> items, ProgressListener listener) throws InterruptedException {
        addAll(items, Runtime.getRuntime().availableProcessors(), listener, DEFAULT_PROGRESS_UPDATE_INTERVAL);
    }

    /**
     * Add multiple items to the index. Reports progress to the passed in implementation of {@link ProgressListener}
     * every progressUpdateInterval elements indexed.
     *
     * @param items the items to add to the index
     * @param numThreads number of threads to use for parallel indexing
     * @param listener listener to report progress to
     * @param progressUpdateInterval after indexing this many items progress will be reported. The last element will always be reported regardless of this setting.
     * @throws InterruptedException thrown when the thread doing the indexing is interrupted
     */
    default void addAll(Collection<TItem> items, int numThreads, ProgressListener listener, int progressUpdateInterval)
            throws InterruptedException {
        ThreadPoolExecutor executorService = new ThreadPoolExecutor(numThreads, numThreads, 60L, TimeUnit.MILLISECONDS,
                new LinkedBlockingQueue<>(),
                new NamedThreadFactory("indexer-%d"));
        executorService.allowCoreThreadTimeOut(true);

        int numItems = items.size();

        AtomicInteger workDone = new AtomicInteger();

        try {
            Queue<TItem> queue = new LinkedBlockingQueue<>(items);
            List<Future<?>> futures = new ArrayList<>();

            for (int threadId = 0; threadId < numThreads; threadId++) {
                futures.add(executorService.submit(() -> {
                    TItem item;
                    while((item = queue.poll()) != null) {
                        add(item);

                        int done = workDone.incrementAndGet();
                        if (done % progressUpdateInterval == 0 || numItems == done) {
                            listener.updateProgress(done, items.size());
                        }
                    }
                }));
            }

            for(Future<?> future : futures) {
                try {
                    future.get();
                } catch (ExecutionException e) {
                    throw new UncategorizedIndexException("An exception was thrown by one of the threads.", e.getCause());
                }
            }

        } finally {
            executorService.shutdown();
        }
    }
    
    /**
     * Add multiple items to the index and, in parallel, searches for similar pairs above a threshold.
     *
     * @param items the items to add to the index
     * @param threshold the threshold used in the search
     * @throws InterruptedException thrown when the thread doing the indexing is interrupted
     */
    //default List<String> addAllJoin(Collection<TItem> items, TDistance threshold) throws InterruptedException {
    default SearchResultHops addAllJoin(Collection<TItem> items, TDistance threshold) throws InterruptedException {
        return addAllJoin(items, threshold, NullProgressListener.INSTANCE);
    }

    /**
     * Add multiple items to the index and, in parallel, searches for similar pairs above a threshold.
     * Reports progress to the passed in implementation of {@link ProgressListener}
     * every {@link Index#DEFAULT_PROGRESS_UPDATE_INTERVAL} elements indexed.
     *
     * @param items the items to add to the index
     * @param threshold the threshold used in the search
     * @param listener listener to report progress to
     * @throws InterruptedException thrown when the thread doing the indexing is interrupted
     */
    //default List<String> addAllJoin(Collection<TItem> items, TDistance threshold, ProgressListener listener) throws InterruptedException {
    default SearchResultHops addAllJoin(Collection<TItem> items, TDistance threshold, ProgressListener listener) throws InterruptedException {
        return addAllJoin(items, threshold, Runtime.getRuntime().availableProcessors(), listener, DEFAULT_PROGRESS_UPDATE_INTERVAL);
    }

    /**
     * Add multiple items to the index and, in parallel, searches for similar pairs above a threshold.
     * Reports progress to the passed in implementation of {@link ProgressListener}
     * every progressUpdateInterval elements indexed.
     *
     * @param items the items to add to the index
     * @param threshold the threshold used in the search
     * @param numThreads number of threads to use for parallel indexing
     * @param listener listener to report progress to
     * @param progressUpdateInterval after indexing this many items progress will be reported. The last element will always be reported regardless of this setting.
     * @throws InterruptedException thrown when the thread doing the indexing is interrupted
     */
    //default List<String> addAllJoin(Collection<TItem> items, TDistance threshold, int numThreads, ProgressListener listener, int progressUpdateInterval)
    default SearchResultHops addAllJoin(Collection<TItem> items, TDistance threshold, int numThreads, ProgressListener listener, int progressUpdateInterval)
            throws InterruptedException {
        ThreadPoolExecutor executorService = new ThreadPoolExecutor(numThreads, numThreads, 60L, TimeUnit.MILLISECONDS,
                new LinkedBlockingQueue<>(),
                new NamedThreadFactory("indexer-%d"));
        executorService.allowCoreThreadTimeOut(true);

        int numItems = items.size();

        AtomicInteger workDone = new AtomicInteger();
        
        AtomicLong totalVisitedVertices = new AtomicLong();
        
        AtomicLong totalHops = new AtomicLong();
        
        //List<String> results = new ArrayList<>();
        ConcurrentLinkedQueue<String> results = new ConcurrentLinkedQueue<>();

        try {
            Queue<TItem> queue = new LinkedBlockingQueue<>(items);
            List<Future<?>> futures = new ArrayList<>();
            
            //numThreads = numThreads/2; // verify
            for (int threadId = 0; threadId < numThreads; threadId++) {
                futures.add(executorService.submit(() -> {
                    TItem item;
                    while((item = queue.poll()) != null) {
                    	
                    	AtomicLong visitedVertices = new AtomicLong();
                    	
                    	AtomicLong hops = new AtomicLong();
                    	
                    	List<String> nodeResults = addJoin(item,threshold, visitedVertices, hops);
                    	
                    	totalVisitedVertices.addAndGet(visitedVertices.get());
                    	
                    	totalHops.addAndGet(hops.get());
                    	
                    	//results.addAll(nodeResults);
                    	synchronized (results) {
                            results.addAll(nodeResults);
                        }
                    	nodeResults = null;
                        int done = workDone.incrementAndGet();
                        if (done % progressUpdateInterval == 0 || numItems == done) {
                            listener.updateProgress(done, items.size());
                        }
                    }
                }));
            }

            for(Future<?> future : futures) {
                try {
                    future.get();
                } catch (ExecutionException e) {
                    throw new UncategorizedIndexException("An exception was thrown by one of the threads.", e.getCause());
                }
            }
            //return results;
            //return new ArrayList<>(results);
            return new SearchResultHops(totalVisitedVertices.get(), totalHops.get(), new ArrayList<>(results));

        } finally {
        	
            executorService.shutdown();
        }
    }
    
    
    /**
     * Add multiple items to the index and, in parallel, searches for similar pairs above a threshold.
     *
     * @param items the items to add to the index
     * @param threshold the threshold used in the search
     * @throws InterruptedException thrown when the thread doing the indexing is interrupted
     */
    default SearchResultTotal addAllJoinTotal(Collection<TItem> items, TDistance threshold) throws InterruptedException {
        return addAllJoinTotal(items, threshold, NullProgressListener.INSTANCE);
    }

    /**
     * Add multiple items to the index and, in parallel, searches for similar pairs above a threshold.
     * Reports progress to the passed in implementation of {@link ProgressListener}
     * every {@link Index#DEFAULT_PROGRESS_UPDATE_INTERVAL} elements indexed.
     *
     * @param items the items to add to the index
     * @param threshold the threshold used in the search
     * @param listener listener to report progress to
     * @throws InterruptedException thrown when the thread doing the indexing is interrupted
     */
    default SearchResultTotal addAllJoinTotal(Collection<TItem> items, TDistance threshold, ProgressListener listener) throws InterruptedException {
        return addAllJoinTotal(items, threshold, Runtime.getRuntime().availableProcessors(), listener, DEFAULT_PROGRESS_UPDATE_INTERVAL);
    }

    /**
     * Add multiple items to the index and, in parallel, searches for similar pairs above a threshold.
     * Reports progress to the passed in implementation of {@link ProgressListener}
     * every progressUpdateInterval elements indexed.
     *
     * @param items the items to add to the index
     * @param threshold the threshold used in the search
     * @param numThreads number of threads to use for parallel indexing
     * @param listener listener to report progress to
     * @param progressUpdateInterval after indexing this many items progress will be reported. The last element will always be reported regardless of this setting.
     * @throws InterruptedException thrown when the thread doing the indexing is interrupted
     */
    default SearchResultTotal addAllJoinTotal(Collection<TItem> items, TDistance threshold, int numThreads, ProgressListener listener, int progressUpdateInterval)
            throws InterruptedException {
        ThreadPoolExecutor executorService = new ThreadPoolExecutor(numThreads, numThreads, 60L, TimeUnit.MILLISECONDS,
                new LinkedBlockingQueue<>(),
                new NamedThreadFactory("indexer-%d"));
        executorService.allowCoreThreadTimeOut(true);

        int numItems = items.size();

        AtomicInteger workDone = new AtomicInteger();
        
        AtomicLong totalVisitedVertices = new AtomicLong();
        
        AtomicLong totalHops = new AtomicLong();
        
        //List<String> results = new ArrayList<>();
        //ConcurrentLinkedQueue<String> results = new ConcurrentLinkedQueue<>();
        
        AtomicLong totalPairs = new AtomicLong();

        try {
            Queue<TItem> queue = new LinkedBlockingQueue<>(items);
            List<Future<?>> futures = new ArrayList<>();
            
            //numThreads = numThreads/2; // verify
            for (int threadId = 0; threadId < numThreads; threadId++) {
                futures.add(executorService.submit(() -> {
                    TItem item;
                    while((item = queue.poll()) != null) {
                    	
                    	AtomicLong visitedVertices = new AtomicLong();
                    	
                    	AtomicLong hops = new AtomicLong();
                    	
                    	List<String> nodeResults = addJoin(item,threshold,visitedVertices,hops);
                    	//results.addAll(nodeResults);
                    	
                    	totalVisitedVertices.addAndGet(visitedVertices.get());
                    	
                    	totalHops.addAndGet(hops.get());
                    	
                    	totalPairs.getAndAdd(nodeResults.size());
                    	nodeResults = null;
                        int done = workDone.incrementAndGet();
                        if (done % progressUpdateInterval == 0 || numItems == done) {
                            listener.updateProgress(done, items.size());
                        }
                    }
                }));
            }

            for(Future<?> future : futures) {
                try {
                    future.get();
                } catch (ExecutionException e) {
                    throw new UncategorizedIndexException("An exception was thrown by one of the threads.", e.getCause());
                }
            }
            //return results;
            //return new ArrayList<>(results);
           // return totalPairs;
            return new SearchResultTotal(totalPairs.get(), 0, totalVisitedVertices.get(), totalHops.get());

        } finally {
        	
            executorService.shutdown();
        }
    }
    
    /**
     * Search multiple items in the the index and, in parallel, for similar pairs above a threshold.
     *
     * @param items the items to add to the index
     * @param threshold the threshold used in the search
     * @throws InterruptedException thrown when the thread doing the indexing is interrupted
     */
    //default List<String> searchAllIndex(Collection<TItem> items, TDistance threshold) throws InterruptedException {
    default SearchResultHops searchAllIndex(Collection<TItem> items, TDistance threshold) throws InterruptedException {
        return searchAllIndex(items, threshold, NullProgressListener.INSTANCE);
    }

    /**
     * Search multiple items in the the index and, in parallel, for similar pairs above a threshold.
     * Reports progress to the passed in implementation of {@link ProgressListener}
     * every {@link Index#DEFAULT_PROGRESS_UPDATE_INTERVAL} elements indexed.
     *
     * @param items the items to add to the index
     * @param threshold the threshold used in the search
     * @param listener listener to report progress to
     * @throws InterruptedException thrown when the thread doing the indexing is interrupted
     */
    //default List<String> searchAllIndex(Collection<TItem> items, TDistance threshold, ProgressListener listener) throws InterruptedException {
    default SearchResultHops searchAllIndex(Collection<TItem> items, TDistance threshold, ProgressListener listener) throws InterruptedException {
        return searchAllIndex(items, threshold, Runtime.getRuntime().availableProcessors(), listener, DEFAULT_PROGRESS_UPDATE_INTERVAL);
    }

    /**
     * Search multiple items in the the index and, in parallel, for similar pairs above a threshold.
     * Reports progress to the passed in implementation of {@link ProgressListener}
     * every progressUpdateInterval elements indexed.
     *
     * @param items the items to add to the index
     * @param threshold the threshold used in the search
     * @param numThreads number of threads to use for parallel indexing
     * @param listener listener to report progress to
     * @param progressUpdateInterval after indexing this many items progress will be reported. The last element will always be reported regardless of this setting.
     * @throws InterruptedException thrown when the thread doing the indexing is interrupted
     */
    //default List<String> searchAllIndex(Collection<TItem> items, TDistance threshold, int numThreads, ProgressListener listener, int progressUpdateInterval)
    default SearchResultHops searchAllIndex(Collection<TItem> items, TDistance threshold, int numThreads, ProgressListener listener, int progressUpdateInterval)
            throws InterruptedException {
        ThreadPoolExecutor executorService = new ThreadPoolExecutor(numThreads, numThreads, 60L, TimeUnit.MILLISECONDS,
                new LinkedBlockingQueue<>(),
                new NamedThreadFactory("indexer-%d"));
        executorService.allowCoreThreadTimeOut(true);

        int numItems = items.size();

        AtomicInteger workDone = new AtomicInteger();
        
        AtomicLong totalVisitedVertices = new AtomicLong();
        
        AtomicLong totalHops = new AtomicLong();
        
        //List<String> results = new ArrayList<>();
        ConcurrentLinkedQueue<String> results = new ConcurrentLinkedQueue<>();

        try {
            Queue<TItem> queue = new LinkedBlockingQueue<>(items);
            List<Future<?>> futures = new ArrayList<>();
            
            //numThreads = 1;//numThreads/2; // verify
            for (int threadId = 0; threadId < numThreads; threadId++) {
                futures.add(executorService.submit(() -> {
                	TItem item;
                    while((item = queue.poll()) != null) {
                    	
                    	AtomicLong visitedVertices = new AtomicLong();
                    	
                    	AtomicLong hops = new AtomicLong();
                    	
                    	List<String> nodeResults = searchParallel(item,threshold, visitedVertices, hops);
                    	
                    	//totalVisitedVertices.addAndGet(visitedVertices.get());
                    	totalVisitedVertices.getAndAdd(visitedVertices.get());
                    	
                    	
                    	//totalHops.addAndGet(hops.get());
                    	totalHops.getAndAdd(hops.get());
                    	
                    	//results.addAll(nodeResults);
                    	synchronized (results) {
                            results.addAll(nodeResults);
                        }
                    	nodeResults = null;
                        int done = workDone.incrementAndGet();
                        if (done % progressUpdateInterval == 0 || numItems == done) {
                            listener.updateProgress(done, items.size());
                        }
                    }
                }));
            }

            for(Future<?> future : futures) {
                try {
                    future.get();
                } catch (ExecutionException e) {
                    throw new UncategorizedIndexException("An exception was thrown by one of the threads.", e.getCause());
                }
            }
            //return results;
            //return new ArrayList<>(results);
            return new SearchResultHops(totalVisitedVertices.get(), totalHops.get(), new ArrayList<>(results));

        } finally {
        	
            executorService.shutdown();
        }
    }
    
    
    /**
     * Search multiple items in the index and, in parallel, retrieve the top-k most similar items.
     *
     * @param items the items to search in the index
     * @param k the number of top similar items to retrieve for each item
     * @throws InterruptedException thrown when the thread performing the search is interrupted
     */
    default SearchResultHops searchAllIndexTopK(Collection<TItem> items, int k, TDistance threshold) throws InterruptedException {
        return searchAllIndexTopK(items, k, threshold, NullProgressListener.INSTANCE);
    }

    /**
     * Search multiple items in the index and, in parallel, retrieve the top-k most similar items.
     * Reports progress to the passed-in implementation of {@link ProgressListener}
     * every {@link Index#DEFAULT_PROGRESS_UPDATE_INTERVAL} elements processed.
     *
     * @param items the items to search in the index
     * @param k the number of top similar items to retrieve for each item
     * @param listener listener to report progress to
     * @throws InterruptedException thrown when the thread performing the search is interrupted
     */
    default SearchResultHops searchAllIndexTopK(Collection<TItem> items, int k, TDistance threshold, ProgressListener listener) throws InterruptedException {
        return searchAllIndexTopK(items, k, threshold, Runtime.getRuntime().availableProcessors(), listener, DEFAULT_PROGRESS_UPDATE_INTERVAL);
    }

    
    /**
     * Search multiple items in the index and, in parallel, retrieve the top-k most similar items.
     * Reports progress to the passed-in implementation of {@link ProgressListener}
     * every progressUpdateInterval elements processed.
     *
     * @param items the items to search in the index
     * @param k the number of top similar items to retrieve for each item
     * @param numThreads number of threads to use for parallel searching
     * @param listener listener to report progress to
     * @param progressUpdateInterval after processing this many items, progress will be reported.
     *                                The last element will always be reported regardless of this setting.
     * @throws InterruptedException thrown when the thread performing the search is interrupted
     */

    default SearchResultHops searchAllIndexTopK(Collection<TItem> items, int k, TDistance threshold, int numThreads, ProgressListener listener, int progressUpdateInterval)
            throws InterruptedException {
        ThreadPoolExecutor executorService = new ThreadPoolExecutor(numThreads, numThreads, 60L, TimeUnit.MILLISECONDS,
                new LinkedBlockingQueue<>(),
                new NamedThreadFactory("indexer-%d"));
        executorService.allowCoreThreadTimeOut(true);

        int numItems = items.size();

        AtomicInteger workDone = new AtomicInteger();
        AtomicLong totalVisitedVertices = new AtomicLong();
        AtomicLong totalHops = new AtomicLong();
        ConcurrentLinkedQueue<String> results = new ConcurrentLinkedQueue<>();

        try {
            Queue<TItem> queue = new LinkedBlockingQueue<>(items);
            List<Future<?>> futures = new ArrayList<>();

            for (int threadId = 0; threadId < numThreads; threadId++) {
                futures.add(executorService.submit(() -> {
                    TItem item;
                    while ((item = queue.poll()) != null) {

                        AtomicLong visitedVertices = new AtomicLong();
                        AtomicLong hops = new AtomicLong();
                        
                        int kLocal = k;

                        List<String> nodeResults = searchParallelTopK(item, kLocal, visitedVertices, hops);
                        
                        totalVisitedVertices.getAndAdd(visitedVertices.get());
                        totalHops.getAndAdd(hops.get());
                        
                        while (!nodeResults.isEmpty()) {
                        	
                            // Take the distance from the last element (less similar)
                            String last = nodeResults.get(0);
                            
                            String[] parts = last.split(",");
                            
                            double distK = Double.parseDouble(parts[2]); 
                            
                            //System.out.println(kLocal + " - " + distK);
                            
                            

                            if (distK < (double) threshold) {
                                break;
                            }
                            
                            //incremento local do visitedVertices, a cada busca
                            AtomicLong v2 = new AtomicLong();
                        	
                            //incremento local do hops, a cada busca
                        	AtomicLong h2 = new AtomicLong();
                            
                            // Heuristic 1: Proportional increment between the similarity of the most distant element and the threshold
							//k = k + (int) Math.round( ((1 - ths[0]) * k) / (1 - simK) );
                            
                            
                           // double simK = 1 - distK;
                            kLocal = kLocal + (int) Math.round(((1 - (double) threshold) * kLocal) / (1 - distK));
                  
                            // Refaz a busca com novo k
                            nodeResults = searchParallelTopK(item, kLocal, v2, h2);
                            
                            totalVisitedVertices.getAndAdd(v2.get());
                            totalHops.getAndAdd(h2.get());
                            
                            //System.out.println(kLocal);
                               
                        }
                        
                        nodeResults.removeIf(result -> {
                            String[] parts = result.split(",");
                            double dist = Double.parseDouble(parts[2]);
                            int id1 = Integer.parseInt(parts[0]);
                            int id2 = Integer.parseInt(parts[1]);
                            return dist < (double) threshold || id1 >= id2;
                        });
                        
                        synchronized (results) {
                            results.addAll(nodeResults);
                        }

                        nodeResults = null;
                        int done = workDone.incrementAndGet();
                        if (done % progressUpdateInterval == 0 || numItems == done) {
                            listener.updateProgress(done, items.size());
                        }
                    }
                }));
            }

            for (Future<?> future : futures) {
                try {
                    future.get();
                } catch (ExecutionException e) {
                    throw new UncategorizedIndexException("An exception was thrown by one of the threads.", e.getCause());
                }
            }

            return new SearchResultHops(totalVisitedVertices.get(), totalHops.get(), new ArrayList<>(results));
        } finally {
            executorService.shutdown();
        }
    }
    
    
    
    /**
     * Search multiple items in the the index and, in parallel, for similar pairs above a threshold.
     *
     * @param items the items to add to the index
     * @param threshold the threshold used in the search
     * @throws InterruptedException thrown when the thread doing the indexing is interrupted
     */
    default SearchResultTotal searchAllIndexTotal(Collection<TItem> items, TDistance threshold) throws InterruptedException {
        return searchAllIndexTotal(items, threshold, NullProgressListener.INSTANCE);
    }

    /**
     * Search multiple items in the the index and, in parallel, for similar pairs above a threshold.
     * Reports progress to the passed in implementation of {@link ProgressListener}
     * every {@link Index#DEFAULT_PROGRESS_UPDATE_INTERVAL} elements indexed.
     *
     * @param items the items to add to the index
     * @param threshold the threshold used in the search
     * @param listener listener to report progress to
     * @throws InterruptedException thrown when the thread doing the indexing is interrupted
     */
    default SearchResultTotal searchAllIndexTotal(Collection<TItem> items, TDistance threshold, ProgressListener listener) throws InterruptedException {
        return searchAllIndexTotal(items, threshold, Runtime.getRuntime().availableProcessors(), listener, DEFAULT_PROGRESS_UPDATE_INTERVAL);
    }
    
    
    /**
     * Search multiple items in the the index and, in parallel, for similar pairs above a threshold.
     * Reports progress to the passed in implementation of {@link ProgressListener}
     * every progressUpdateInterval elements indexed.
     *
     * @param items the items to add to the index
     * @param threshold the threshold used in the search
     * @param numThreads number of threads to use for parallel indexing
     * @param listener listener to report progress to
     * @param progressUpdateInterval after indexing this many items progress will be reported. The last element will always be reported regardless of this setting.
     * @throws InterruptedException thrown when the thread doing the indexing is interrupted
     */
    default SearchResultTotal searchAllIndexTotal(Collection<TItem> items, TDistance threshold, int numThreads, ProgressListener listener, int progressUpdateInterval)
            throws InterruptedException {
        ThreadPoolExecutor executorService = new ThreadPoolExecutor(numThreads, numThreads, 60L, TimeUnit.MILLISECONDS,
                new LinkedBlockingQueue<>(),
                new NamedThreadFactory("indexer-%d"));
        executorService.allowCoreThreadTimeOut(true);

        int numItems = items.size();

        AtomicInteger workDone = new AtomicInteger();
        
        AtomicLong totalVisitedVertices = new AtomicLong();
        
        AtomicLong totalHops = new AtomicLong();
        
        //List<String> results = new ArrayList<>();
        //ConcurrentLinkedQueue<String> results = new ConcurrentLinkedQueue<>();
        
        AtomicLong totalPairs = new AtomicLong(0);

        try {
            Queue<TItem> queue = new LinkedBlockingQueue<>(items);
            List<Future<?>> futures = new ArrayList<>();
            
            //numThreads = numThreads/2; // verify
            for (int threadId = 0; threadId < numThreads; threadId++) {
                futures.add(executorService.submit(() -> {
                	TItem item;
                    while((item = queue.poll()) != null) {
                    	
                    	AtomicLong visitedVertices = new AtomicLong();
                    	
                    	AtomicLong hops = new AtomicLong();
                    	
                    	List<String> nodeResults = searchParallel(item,threshold,visitedVertices,hops);
                    	
                    	totalVisitedVertices.addAndGet(visitedVertices.get());
                    	
                    	totalHops.addAndGet(hops.get());
                    	
                    	//totalPairs
                    	totalPairs.getAndAdd(nodeResults.size());
                    	
                    	//results.addAll(nodeResults);
                    	nodeResults = null;
                        int done = workDone.incrementAndGet();
                        if (done % progressUpdateInterval == 0 || numItems == done) {
                            listener.updateProgress(done, items.size());
                        }
                    }
                }));
            }

            for(Future<?> future : futures) {
                try {
                    future.get();
                } catch (ExecutionException e) {
                    throw new UncategorizedIndexException("An exception was thrown by one of the threads.", e.getCause());
                }
            }
            //return results;
            //return new ArrayList<>(results);
            
            //return totalPairs
            //return totalPairs;
            return new SearchResultTotal(totalPairs.get(), 0, totalVisitedVertices.get(), totalHops.get());

        } finally {
        	
            executorService.shutdown();
        }
    }
    
    
    /**
     * Search multiple items in the the index and, in parallel, for similar pairs above a threshold.
     *
     * @param items the items to add to the index
     * @param threshold the threshold used in the search
     * @throws InterruptedException thrown when the thread doing the indexing is interrupted
     */
    default SearchResultTotal searchAllIndexTopKTotal(Collection<TItem> items, int k, TDistance threshold) throws InterruptedException {
        return searchAllIndexTopKTotal(items, k, threshold, NullProgressListener.INSTANCE);
    }

    /**
     * Search multiple items in the the index and, in parallel, for similar pairs above a threshold.
     * Reports progress to the passed in implementation of {@link ProgressListener}
     * every {@link Index#DEFAULT_PROGRESS_UPDATE_INTERVAL} elements indexed.
     *
     * @param items the items to add to the index
     * @param threshold the threshold used in the search
     * @param listener listener to report progress to
     * @throws InterruptedException thrown when the thread doing the indexing is interrupted
     */
    default SearchResultTotal searchAllIndexTopKTotal(Collection<TItem> items, int k, TDistance threshold, ProgressListener listener) throws InterruptedException {
        return searchAllIndexTopKTotal(items, k, threshold, Runtime.getRuntime().availableProcessors(), listener, DEFAULT_PROGRESS_UPDATE_INTERVAL);
    }
    
    
    /**
     * Search multiple items in the the index and, in parallel, for similar pairs above a threshold.
     * Reports progress to the passed in implementation of {@link ProgressListener}
     * every progressUpdateInterval elements indexed.
     *
     * @param items the items to add to the index
     * @param threshold the threshold used in the search
     * @param numThreads number of threads to use for parallel indexing
     * @param listener listener to report progress to
     * @param progressUpdateInterval after indexing this many items progress will be reported. The last element will always be reported regardless of this setting.
     * @throws InterruptedException thrown when the thread doing the indexing is interrupted
     */
    default SearchResultTotal searchAllIndexTopKTotal(Collection<TItem> items, int k, TDistance threshold, int numThreads, ProgressListener listener, int progressUpdateInterval)
            throws InterruptedException {
        ThreadPoolExecutor executorService = new ThreadPoolExecutor(numThreads, numThreads, 60L, TimeUnit.MILLISECONDS,
                new LinkedBlockingQueue<>(),
                new NamedThreadFactory("indexer-%d"));
        executorService.allowCoreThreadTimeOut(true);

        int numItems = items.size();

        AtomicInteger workDone = new AtomicInteger();
        
        AtomicLong totalVisitedVertices = new AtomicLong();
        
        AtomicLong totalHops = new AtomicLong();
        
        AtomicLong totalPairs = new AtomicLong(0);

        try {
            Queue<TItem> queue = new LinkedBlockingQueue<>(items);
            List<Future<?>> futures = new ArrayList<>();
            
            for (int threadId = 0; threadId < numThreads; threadId++) {
                futures.add(executorService.submit(() -> {
                	TItem item;
                    while((item = queue.poll()) != null) {
                    	
                    	AtomicLong visitedVertices = new AtomicLong();
                    	
                    	AtomicLong hops = new AtomicLong();
                    	
                    	int kLocal = k;
                    	
                    	List<String> nodeResults = searchParallelTopK(item, kLocal, visitedVertices, hops);
                    	
                        totalVisitedVertices.getAndAdd(visitedVertices.get());
                        totalHops.getAndAdd(hops.get());
                        
                        //incremento local do visitedVertices, a cada busca
                        AtomicLong v2 = new AtomicLong();
                    	
                        //incremento local do hops, a cada busca
                    	AtomicLong h2 = new AtomicLong();
                    	
                    	while (!nodeResults.isEmpty()) {
                        	
                            // Take the distance from the last element (less similar)
                            String last = nodeResults.get(0);
                            
                            String[] parts = last.split(",");
                            
                            double distK = Double.parseDouble(parts[2]); 
                            
                            //System.out.println(kLocal + " - " + distK);
                            

                            if (distK < (double) threshold) {
                                break;
                            }
                            
                            nodeResults = null;
                            
                            // Heuristic 1: Proportional increment between the similarity of the most distant element and the threshold
							//k = k + (int) Math.round( ((1 - ths[0]) * k) / (1 - simK) );
                            
                            
                           // double simK = 1 - distK;
                            kLocal = kLocal + (int) Math.round(((1 - (double) threshold) * kLocal) / (1 - distK));
                  
                            // Refaz a busca com novo k
                            nodeResults = searchParallelTopK(item, kLocal, v2, h2);
                            
                            totalVisitedVertices.getAndAdd(v2.get());
                            totalHops.getAndAdd(h2.get());
                            
                            v2.set(0);
                            h2.set(0);
                            
                            //System.out.println(kLocal);
                               
                        }                    	
                    	
                    	nodeResults.removeIf(result -> {
                            String[] parts = result.split(",");
                            double dist = Double.parseDouble(parts[2]);
                            int id1 = Integer.parseInt(parts[0]);
                            int id2 = Integer.parseInt(parts[1]);
                            return dist < (double) threshold || id1 >= id2;
                        });
                    	
                    	
                    	
                    	//totalPairs
                    	totalPairs.getAndAdd(nodeResults.size());
                    	
                    	//results.addAll(nodeResults);
                    	nodeResults = null;
                        int done = workDone.incrementAndGet();
                        if (done % progressUpdateInterval == 0 || numItems == done) {
                            listener.updateProgress(done, items.size());
                        }
                    }
                }));
            }

            for(Future<?> future : futures) {
                try {
                    future.get();
                } catch (ExecutionException e) {
                    throw new UncategorizedIndexException("An exception was thrown by one of the threads.", e.getCause());
                }
            }
            //return results;
            //return new ArrayList<>(results);
            
            //return totalPairs
            //return totalPairs;
            return new SearchResultTotal(totalPairs.get(), 0, totalVisitedVertices.get(), totalHops.get());

        } finally {
        	
            executorService.shutdown();
        }
    }
    

    /**
     * Returns the size of the index.
     *
     * @return size of the index
     */
    int size();

    /**
     * Returns an item by its identifier.
     *
     * @param id unique identifier or the item to return
     * @return an item
     */
    Optional<TItem> get(TId id);

    /**
     * Returns all items in the index.
     *
     * @return all items in the index
     */
    Collection<TItem> items();

    /**
     * Find the items closest to the passed in vector.
     *
     * @param vector the vector
     * @param k number of items to return
     * @return the items closest to the passed in vector
     */
    List<SearchResult<TItem, TDistance>> findNearest(TVector vector, int k, AtomicLong visitedVertices, AtomicLong hops);
    
    /**
     * Find the items closest to the passed in vector.
     *
     * @param vector the vector
     * @param k number of items to return
     * @return the items closest to the passed in vector
     */
    SearchResultTotal findNearestTotal(TVector vector, int k, AtomicLong visitedVertices, AtomicLong hops);

    /**
     * Find the items closest to the item identified by the passed in id. If the id does not match an item an empty
     * list is returned. the element itself is not included in the response.
     *
     * @param id id of the item to find the neighbors of
     * @param k number of items to return
     * @return the items closest to the item
     */
    //default List<SearchResult<TItem, TDistance>> findNeighbors(TId id, int k) {
    //   return get(id).map(item -> findNearest(item.vector(), k + 1).stream()
    //            .filter(result -> !result.item().id().equals(id))
    //            .limit(k)
    //            .collect(Collectors.toList()))
    //            .orElse(Collections.emptyList());
    //}
    default List<SearchResult<TItem, TDistance>> findNeighbors(TId id, int k) {
    	AtomicLong visitedVertices = new AtomicLong(0);
    	AtomicLong hops = new AtomicLong(0);
        return get(id).map(item -> findNearest(item.vector(), k + 1, visitedVertices, hops).stream()
                 .filter(result -> !result.item().id().equals(id))
                 .limit(k)
                 .collect(Collectors.toList()))
                 .orElse(Collections.emptyList());
     }

    /**
     * Saves the index to an OutputStream. Saving is not thread safe and you should not modify the index while saving.
     *
     * @param out the output stream to write the index to
     * @throws IOException in case of I/O exception
     */
    void save(OutputStream out) throws IOException;

    /**
     * Saves the index to a file. Saving is not thread safe and you should not modify the index while saving.
     *
     * @param file file to write the index to
     * @throws IOException in case of I/O exception
     */
    default void save(File file) throws IOException {
        save(new FileOutputStream(file));
    }

    /**
     * Saves the index to a path. Saving is not thread safe and you should not modify the index while saving.
     *
     * @param path file to write the index to
     * @throws IOException in case of I/O exception
     */
    default void save(Path path) throws IOException {
        save(Files.newOutputStream(path));
    }

}