import java.util.*;
import java.io.*;
import org.javatuples.Pair;
import org.javatuples.Triplet;
import org.javatuples.Quartet;

/**
 * The algorithms used to find Simpson's Paradoxes in multi-dimensional data sets.
 *
 * @Jay Xu
 * @July 2024
 */
public class Redundancy
{
    public static final int size = 8; // number of attributes
    public static final String filename = "adult.csv";
    public static final boolean strong = true;
    public static long time = System.currentTimeMillis();
    public static int SPCounter;
    public static int nonRedundantCounter;
    /**
     * parses the csv
     * 
     */
    public static List<String[]> getLines(String fileName) {
        String line;
        List<String[]> records = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(fileName))) {
            while ((line = br.readLine()) != null) {
                String[] fields = line.split(",", 0);
                boolean flag = false;
                for(String s : fields) {
                    if(s.contains("?")) {
                        flag = true;
                        break;
                    }
                }
                if(flag){
                    continue;
                }
                records.add(fields);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return records;
    }
    public static void createIndex(List<String[]> records, HashMap<String, Integer> indexHolder, int[] lengths, List<List<Integer>> enumRecs, List<Integer> binaryLabels) {
        boolean firstLine = true;
        for (String[] record : records) {
            if(firstLine){ 
                firstLine = false;
                continue;
            }
            int count = 0;
            List<Integer> enumRec = new ArrayList<>();
            for(String s: record){
                if(count == size){
                    binaryLabels.add(Integer.parseInt(s));
                    break;
                }
                // create a distinct key that corresponds with an index
                String key = count + " " + s; 
                if(!indexHolder.containsKey(key)){
                    // update the amount of elements for each attribute.
                   lengths[count]++; 
                   int i = lengths[count];
                indexHolder.put(key, i); 
                
                }
                enumRec.add(indexHolder.get(key));
                count++;
            }
            enumRecs.add(enumRec);

        }
    } 

    public static int[] getMover(int[] lengths){ 
        int[] bitSize = new int[size];
        bitSize[0] = 0;
        for(int i = 1; i < size; i++){
            bitSize[i] = bitSize[i-1] + (int)(Math.log(lengths[i-1])/Math.log(2)) + 1; // how many bits to move each attribute over
            
        }
        
        return bitSize;
        
    }
    public static void aggregate(int[] move,HashMap<String, Integer> indexHolder, HashMap<Long, Long> map, List<List<Integer>> enumRecs, List<Integer> binaryLabels,  HashMap<Long, List<Long>> covers) {
        HashMap<List<Integer>, Long> agg = new HashMap<>();
        int count = 0;
        for(List<Integer> listKey : enumRecs) {
            if(!agg.containsKey(listKey)) {
                agg.put(listKey, (long)0);
            }
            agg.put(listKey, agg.get(listKey) + ((long)1 << 32) +(long)binaryLabels.get(count));
            count++;
        }
        for(List<Integer> listKey: agg.keySet()) {
            long base = 0;
            for(int j = 0; j < size; j++) {
                base += (long)listKey.get(j) << move[j];
            }
            for(int i = 0; i < (1 << size); i++) {
                long l = 0;
                for(int j = 0; j < size; j++){
                    long ind = (long) 1<< j;
                    if((ind & i) == ind) {
                        l += (long)listKey.get(j) << move[j];
                    }
                }
                if(!map.containsKey(l)){
                    map.put(l, (long)0);
                    covers.put(l, new ArrayList<>());
                }
                map.put(l, map.get(l) + agg.get(listKey));
                covers.get(l).add(base);
            }

        }



    }
    public static double getRate(HashMap<Long, Long> map, Long key){ 
        double rate = 0;
        if(map.containsKey(key)){
            Long l1 = map.get(key);
            long l2 =  l1.longValue();
            double inc = l2 % ((long)1 << 32);
            double tot = l2 >> 32;
            rate = inc / tot;  
        }
        return rate;
    }

    public static HashMap<List<Long>, List<Long>> groupCovers(HashMap<Long, List<Long>> covers) {
        HashMap<List<Long>, List<Long>> gc = new HashMap<>();
        for(Long pop: covers.keySet()) {
            List<Long> newKey = covers.get(pop);
            if(!gc.containsKey(newKey)) {
                gc.put(newKey, new ArrayList<>());
            }
            gc.get(newKey).add(pop);
        }
        return gc;

    }
    public static boolean checkSP(Long c1, Long c2, int j, HashMap<Long, Long> stats, int[] move, int[] lengths) {
        if(!stats.containsKey(c1) || !stats.containsKey(c2)) {
            return false;
        }

        if((stats.get(c1) >> 32) + (stats.get(c2) >> 32) <= 4* lengths[j]) {
            return false;
        }
        
        double c1r = getRate(stats, c1);
        double c2r = getRate(stats, c2);

        if(c1r > c2r) {
            for(long k = 1; k <= lengths[j]; k++){
                long c1x = c1 + (k << move[j]);
                long c2x = c2 + (k << move[j]);
                if(!stats.containsKey(c1x) || !stats.containsKey(c2x)) {
                    return false;
                }
                double c1xr = getRate(stats, c1x);
                double c2xr = getRate(stats, c2x);

                if(c1xr > c2xr){
                    return false;
                }

            }
        } else if (c1r < c2r) {
            for(long k = 1; k <= lengths[j]; k++){
                long c1x = c1 + (k << move[j]);
                long c2x = c2 + (k << move[j]);
                if(!stats.containsKey(c1x) || !stats.containsKey(c2x)) {
                    return false;
                }
                double c1xr = getRate(stats, c1x);
                double c2xr = getRate(stats, c2x);

                if(c1xr < c2xr){
                    return false;
                }

            }
        } else {
            int c1gc2 = 0, c2gc1 = 0;
            for(long k = 1; k <= lengths[j]; k++){
                long c1x = c1 + (k << move[j]);
                long c2x = c2 + (k << move[j]);
                if(!stats.containsKey(c1x) || !stats.containsKey(c2x)) {
                    return false;
                }
                double c1xr = getRate(stats, c1x);
                double c2xr = getRate(stats, c2x);

                if(c1xr < c2xr){
                    c2gc1++;

                } else if(c1xr > c2xr) {
                    c1gc2++;
                }

                

            }

            return (c1gc2 == 0 || c2gc1 == 0) && (c1gc2 != 0 || c2gc1 != 0);

        }
        return true;
    }

    public static Pair<List<Integer>, Long> getLeastStars(List<Long> list, int[] mover) {
        int least = size + 1;
        long key = 0;
        List<Integer> zeros = new ArrayList<>();
        for(Long l: list){
            int numZeros = 0;
            List<Integer> nzeros = new ArrayList<>();
            for(int j = 0; j < size - 1; j++) {
                
                if(((l >> mover[j]) & (((long) 1 << (mover[j+1] - mover[j]) ) - 1)) == 0) {
                    numZeros++;
                    nzeros.add(j);
                }


            }
            if(numZeros < least) {
                key = l;
                least = numZeros;
                zeros = nzeros;
            }
        }
        return new Pair<>(zeros, key);
    }
    public static List<Triplet<Long, Long, Integer>> findSP(HashMap<Long,Long> aggregations, HashMap<List<Long>, List<Long>> groupedCovers, int[] lengths, int[] mover ){
        List<Triplet<Long, Long, Integer>> allSps = new ArrayList<>(); // c1, c2, j
        for(List<Long> k : groupedCovers.keySet()) { // k is the covering, v is the list of populations that have that covering
            List<Long> v = groupedCovers.get(k);
            Pair<List<Integer>, Long> info = getLeastStars(v, mover);
            long l1 = info.getValue1();
            List<Integer> zeros = info.getValue0();
            List<Quartet<Integer, Integer, Long, Long>> coverSPs = new ArrayList<>(); // i, j, i1, i2
            if(zeros.size() < 2){
                continue;
            }
            for(int i : zeros){
                for(long i1 = 1; i1 < lengths[i]; i1++){
                    long c1 = l1 + (i1 << mover[i]);
                    
                    if(!aggregations.containsKey(c1)){
                        continue;
                    }
                    for(long i2 = i1 + 1; i2 <= lengths[i]; i2++){
                        long c2 = l1 + (i2 << mover[i]);
                        if(!aggregations.containsKey(c2)){
                            continue;
                        }
                        
                        for(int j :zeros){
                            if(i == j) {
                                continue;
                            }

                            boolean isSP = checkSP(c1, c2, j, aggregations, mover, lengths);
                            if(isSP) {
                                coverSPs.add(new Quartet<>(i, j, i1, i2));
                            }
                        }
                    }

                }

            }
            for(Quartet<Integer, Integer, Long, Long> q : coverSPs){
                nonRedundantCounter++;
                int i = q.getValue0();
                int j = q.getValue1();
                long i1 = q.getValue2();
                long i2 = q.getValue3();
                for(long pop:v){ // for all the populations that have this covering, call it an sp
                    Triplet<Long, Long, Integer> t = new Triplet<>(pop + (i1 << mover[i]), pop + (i2 << mover[i]), j);
                    allSps.add(t);
                    SPCounter++;
                }

            }
        }
        return allSps;
    }


    /**
     * Example of Algorithm 1
     */
    public static void main(String[] args){
        HashMap<String, Integer> indexHolder = new HashMap<>(); // Map to hold indexes
       
        int[] lengths = new int[size];// the 
        List<String[]> lines = getLines(filename);
        List<List<Integer>> enumRecs = new ArrayList<>();
        List<Integer> bLabs = new ArrayList<>();
        createIndex(lines, indexHolder, lengths, enumRecs, bLabs);
        int[] move = getMover(lengths);
        HashMap<Long,Long> aggregations = new HashMap<>();
        
        HashMap<Long, List<Long>> covers = new HashMap<>();
        
        //System.out.println(System.currentTimeMillis() - time);
        aggregate(move, indexHolder, aggregations, enumRecs, bLabs, covers);
        HashMap<List<Long>, List<Long>> groupedCovers = groupCovers(covers);
        System.out.println(covers.size());
        System.out.println(groupedCovers.size());
        System.out.println(System.currentTimeMillis() - time);
        SPCounter = 0;
        nonRedundantCounter = 0;
        List<Triplet<Long, Long, Integer>> SPs = findSP(aggregations, groupedCovers, lengths, move);

        
        System.out.println(SPCounter);
        System.out.println(nonRedundantCounter);
        System.out.println(System.currentTimeMillis() - time);
    }
  
}