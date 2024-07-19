import java.util.*;
import java.io.*;
import org.javatuples.Pair;
import org.javatuples.Triplet;
import org.javatuples.Quartet;

public class Redundancy2 {
    public static final int size = 8; // number of attributes (adult : 8) (loan : 8) (mushroom : 12)
    public static final String filename = "loan.csv";
    public static final boolean strong = true;
    public static long time = System.currentTimeMillis();
    public static int SPCounter;
    public static int nonRedundantCounter;
    public static int smth = 0;
    /**
     * parses the csv
     * @param fileName
     * @return a list of the records
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
                if (fields[size].equals("e")) {
                    fields[size] = "1";
                } else if (fields[size].equals("p")) {
                    fields[size] = "0";
                } else {
                    smth = 0;
                }
                records.add(fields);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return records;
    }

    /**
     * Creates the indexes and enumerates the records
     * @param records string list of records
     * @param indexHolder index holder that maps a value to an index
     * @param lengths the length of each variable
     * @param enumRecs the enumerated list of records
     * @param binaryLabels the list of binary labels
     */
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

    /**
     * Calculates the amount of bits we need to shift for each index
     * @param lengths
     * @return int[] mover
     */
    public static int[] getMover(int[] lengths){
        int[] bitSize = new int[size];
        bitSize[0] = 0;
        for(int i = 1; i < size; i++){
            bitSize[i] = bitSize[i-1] + (int)(Math.log(lengths[i-1])/Math.log(2)) + 1; // how many bits to move each attribute over

        }

        return bitSize;

    }

    /**
     * Gets the stats and covers
     * @param move
     * @param indexHolder the index holder
     * @param map the map we will put the stats in
     * @param enumRecs the records in enumerated form
     * @param binaryLabels the binary labels corresponding with records
     * @param covers the covers of each population
     */
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

    /**
     * gets the rate for a population
     * @param map
     * @param key the population
     * @return the rate of the population
     */
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

    /**
     * Groups covers together so that the covering is  now the key and the value is list of populations with the same cover
     * @param covers
     * @return grouped covering
     */
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

    /**
     * Checks whether the two siblings and separator are a SP
     * @param c1 sibling 1
     * @param c2 sibling 2
     * @param j separator
     * @param stats the stats of the populations
     * @param move
     * @param lengths
     * @return boolean whether it is an SP
     */
    public static Pair<Boolean, Set<List<Long>>> checkSP(Long c1, Long c2, int j, HashMap<Long, Long> stats, int[] move, int[] lengths, HashMap<Long, List<Long>> covers) {
        Set<List<Long>> separatorsInvolved = new HashSet<>();
        if(!stats.containsKey(c1) || !stats.containsKey(c2)) {
            return new Pair<Boolean, Set<List<Long>>>(false, separatorsInvolved);
        }

        separatorsInvolved.add(covers.get(c1));
        separatorsInvolved.add(covers.get(c2));

        double c1r = getRate(stats, c1);
        double c2r = getRate(stats, c2);

        if(c1r > c2r) {
            for(long k = 1; k <= lengths[j]; k++){
                long c1x = c1 + (k << move[j]);
                long c2x = c2 + (k << move[j]);
                if(!stats.containsKey(c1x) || !stats.containsKey(c2x)) {
                    return new Pair<Boolean, Set<List<Long>>>(false, separatorsInvolved);
                }
                double c1xr = getRate(stats, c1x);
                double c2xr = getRate(stats, c2x);

                if(c1xr > c2xr){
                    return new Pair<Boolean, Set<List<Long>>>(false, separatorsInvolved);
                }
                separatorsInvolved.add(covers.get(c1x));
                separatorsInvolved.add(covers.get(c2x));


            }
        } else if (c1r < c2r) {
            for(long k = 1; k <= lengths[j]; k++){
                long c1x = c1 + (k << move[j]);
                long c2x = c2 + (k << move[j]);
                if(!stats.containsKey(c1x) || !stats.containsKey(c2x)) {
                    return new Pair<Boolean, Set<List<Long>>>(false, separatorsInvolved);
                }
                double c1xr = getRate(stats, c1x);
                double c2xr = getRate(stats, c2x);

                if(c1xr < c2xr){
                    return new Pair<Boolean, Set<List<Long>>>(false, separatorsInvolved);
                }

                separatorsInvolved.add(covers.get(c1x));
                separatorsInvolved.add(covers.get(c2x));

            }
        } else {
            int c1gc2 = 0, c2gc1 = 0;
            for(long k = 1; k <= lengths[j]; k++){
                long c1x = c1 + (k << move[j]);
                long c2x = c2 + (k << move[j]);
                if(!stats.containsKey(c1x) || !stats.containsKey(c2x)) {
                    return new Pair<Boolean, Set<List<Long>>>(false, separatorsInvolved);
                }
                double c1xr = getRate(stats, c1x);
                double c2xr = getRate(stats, c2x);
                if(c1xr < c2xr){
                    c2gc1++;
                } else if(c1xr > c2xr) {
                    c1gc2++;
                } else {
                    return new Pair<Boolean, Set<List<Long>>>(false, separatorsInvolved);
                }

                separatorsInvolved.add(covers.get(c1x));
                separatorsInvolved.add(covers.get(c2x));
            }
            // return (c1gc2 == 0 || c2gc1 == 0) && (c1gc2 != 0 || c2gc1 != 0);
            return new Pair<Boolean, Set<List<Long>>>((c1gc2 == 0 || c2gc1 == 0), separatorsInvolved);

        }
        // return true;
        return new Pair<Boolean, Set<List<Long>>>(true, separatorsInvolved);
    }

    /**
     * Gets the population with the least stars
     * @param list of populations
     * @param mover
     * @return
     */
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
            if((l >> mover[size - 1]) == 0) {
                numZeros++;
                nzeros.add(size - 1);
            }
            if(numZeros < least) {
                key = l;
                least = numZeros;
                zeros = nzeros;
            }
        }
        return new Pair<>(zeros, key);
    }

    public static Long[] getIndexes(Long key, int[] mover){
        Long[] ret = new Long[size];
        for(int j = 0; j < size - 1; j++) {
            long l = ((key >> mover[j]) & (((long) 1 << (mover[j+1] - mover[j]) ) - 1));
            ret[j] = l;
        }
        ret[size - 1] = key >> mover[size-1];
        return ret;
    }

    /**
     * Finds the list of SP
     * @param aggregations the statistics of each population
     * @param groupedCovers the grouped covers
     * @param lengths
     * @param mover
     * @return a list of all SP
     */
    public static Pair<Set<Triplet<Long, Long, Integer>>, HashMap<Set<List<Long>>, Set<Triplet<Long, Long, Integer>>>> findSP(HashMap<Long,Long> aggregations, HashMap<List<Long>, List<Long>> groupedCovers, int[] lengths, int[] mover, HashMap<Long, List<Long>> covers ){
        Set<Triplet<Long, Long, Integer>> allSps = new HashSet<>(); // c1, c2, j
        HashMap<Set<List<Long>>, Set<Triplet<Long, Long, Integer>>> infos = new HashMap<>();
        HashMap<Set<List<Long>>, Set<Triplet<Long, Long, Integer>>> infosAux = new HashMap<>();

        for(List<Long> k : groupedCovers.keySet()) { // k is the covering, v is the list of populations that have that covering
            List<Long> v = groupedCovers.get(k);
            Pair<List<Integer>, Long> stars = getLeastStars(v, mover); // use the population with least stars
            long l1 = stars.getValue1();
            List<Integer> zeros = stars.getValue0();
            // List<Quartet<Integer, Integer, Long, Long>> coverSPs = new ArrayList<>(); // i, j, i1, i2
            // List<Set<List<Long>>> siblingCovers = new ArrayList<>();
            if(zeros.size() < 2){
                continue;
            }
            for(int i : zeros){ // for the stars
                for(long i1 = 1; i1 < lengths[i]; i1++){ // for the indexes in that star to make siblings
                    long c1 = l1 + (i1 << mover[i]);

                    if(!aggregations.containsKey(c1)){
                        continue;
                    }
                    for(long i2 = i1 + 1; i2 <= lengths[i]; i2++){ // sibling 2
                        long c2 = l1 + (i2 << mover[i]);
                        if(!aggregations.containsKey(c2)){
                            continue;
                        }

                        for(int j : zeros){ // for seperators
                            if(i == j) {
                                continue;
                            }

                            Triplet<Long, Long, Integer> p = new Triplet<>(c1, c2, j);
                            if(allSps.contains(p)){
                                continue;
                            }

                            Pair<Boolean, Set<List<Long>>> res = checkSP(c1, c2, j, aggregations, mover, lengths , covers); // check SP
                            boolean isSp = res.getValue0();
                            Set<List<Long>> info = res.getValue1();
                            if(isSp && !infosAux.containsKey(info)) {
                                allSps.add(p);
                                infosAux.put(info, new HashSet<>());
                                infos.put(info, new HashSet<>());
                                infosAux.get(info).add(p);
                                infos.get(info).add(p);

                                // coverSPs.add(new Quartet<>(i, j, i1, i2)); // if it is an SP, add it to this cover
                                // siblingCovers.add(seperatorsInvolved);
                                // do the I stuff
                                // get all the populations with the same cover

                                List<Long> c1cov = groupedCovers.get(covers.get(c1));
                                List<Long> c2cov = groupedCovers.get(covers.get(c2));
                                for (Long t1 : c1cov) {
                                        Long[] t1Ind = getIndexes(t1, mover);

                                        for (Long t2 : c2cov) {
                                            Long[] t2Ind = getIndexes(t2, mover);
                                            int diff = 0;
                                            for (int ind = 0; ind < size; ind++) {
                                                if (t1Ind[ind] != t2Ind[ind]) {
                                                    diff++;
                                                }
                                            }
                                            if (diff == 1) { // sibings
                                                long small = t1;
                                                long big = t2;
                                                if(t1 > t2){
                                                    small = t2;
                                                    big = t1;
                                                }
                                                Triplet<Long, Long, Integer> p_prime = new Triplet<>(small, big, j);
                                                allSps.add(p_prime);
                                                infos.get(info).add(p_prime);
                                                infosAux.get(info).add(p_prime);
                                            }
                                        }
                                    }
                            } else if (isSp && infosAux.containsKey(info)) {
                                for (Triplet<Long, Long, Integer> p_pprime : infosAux.get(info)) {
                                    long e1 = p_pprime.getValue0();
                                    long e2 = p_pprime.getValue1();
                                    Triplet<Long, Long, Integer> p_hat = new Triplet<>(e1, e2, j);
                                    allSps.add(p_hat);
                                    infos.get(info).add(p_hat);
                                }
                            } else {
                                continue;
                            }
                                // if(!infosAux.containsKey(seperatorsInvolved)) {
                                    // for (Long t1 : c1cov) {
                                    //     Long[] t1Ind = getIndexes(t1, mover);

                                    //     for (Long t2 : c2cov) {
                                    //         Long[] t2Ind = getIndexes(t2, mover);
                                    //         int diff = 0;
                                    //         for (int ind = 0; ind < size; ind++) {
                                    //             if (t1Ind[ind] != t2Ind[ind]) {
                                    //                 diff++;
                                    //             }
                                    //         }
                                    //         if (diff == 1) { // sibings
                                    //             if (!infosAux.containsKey(seperatorsInvolved)) {
                                    //                 infosAux.put(seperatorsInvolved, new HashSet<>());
                                    //                 infos.put(seperatorsInvolved, new HashSet<>());

                                    //             }
                                    //             long small = t1;
                                    //             long big = t2;
                                    //             if(t1 > t2){
                                    //                 small = t2;
                                    //                 big = t1;
                                    //             }
                                    //             Triplet<Long, Long, Integer> p = new Triplet<>(small, big, j);
                                    //             allSps.add(p);
                                    //             infos.get(seperatorsInvolved).add(p);
                                    //             infosAux.get(seperatorsInvolved).add(p);
                                    //         }
                                    //     }
                                    // }
                                // } else {
                                //     for (Triplet<Long, Long, Integer> p : infosAux.get(seperatorsInvolved)) {
                                //         Triplet<Long, Long, Integer> pprime = new Triplet<>(p.getValue0(), p.getValue1(), j);
                                //         allSps.add(pprime);
                                //         infos.get(seperatorsInvolved).add(pprime);
                                //         smth++;
                                //         //System.out.println("This actually triggers");
                                //     }
                                // }
                            // }
                        }
                    }

                }

            }
        }
        return new Pair<Set<Triplet<Long, Long, Integer>>, HashMap<Set<List<Long>>, Set<Triplet<Long, Long, Integer>>>>(allSps, infos);
    }


    /**
     * Example of Algorithm 1
     */
    public static void main(String[] args){
        HashMap<String, Integer> indexHolder = new HashMap<>(); // Map to hold indexes

        int[] lengths = new int[size];//
        System.out.println("Dataset: " + filename);
        List<String[]> lines = getLines(filename);
        System.out.println("Number of Records: " + lines.size());
        List<List<Integer>> enumRecs = new ArrayList<>();
        List<Integer> bLabs = new ArrayList<>();
        createIndex(lines, indexHolder, lengths, enumRecs, bLabs);
        int[] move = getMover(lengths);
        HashMap<Long,Long> aggregations = new HashMap<>();

        HashMap<Long, List<Long>> covers = new HashMap<>();

        // System.out.println(System.currentTimeMillis() - time);
        aggregate(move, indexHolder, aggregations, enumRecs, bLabs, covers);
        HashMap<List<Long>, List<Long>> groupedCovers = groupCovers(covers);
        System.out.println("Total number of non-empty populations: " + covers.size());
        System.out.println("Total number of coverage equivalent subsets: " + groupedCovers.size());
        // System.out.println(System.currentTimeMillis() - time);
        SPCounter = 0;
        nonRedundantCounter = 0;
        // Set<Triplet<Long, Long, Integer>> SPs = findSP(aggregations, groupedCovers, lengths, move, covers);
        Pair<Set<Triplet<Long, Long, Integer>>, HashMap<Set<List<Long>>, Set<Triplet<Long, Long, Integer>>>> tup = findSP(aggregations, groupedCovers, lengths, move, covers);
        Set<Triplet<Long, Long, Integer>> SPs = tup.getValue0();
        HashMap<Set<List<Long>>, Set<Triplet<Long, Long, Integer>>> infos = tup.getValue1();


        System.out.println("Total number of SP: " + SPs.size());
        System.out.println("Total number of redundant relations: " + infos.size());
        // System.out.println(nonRedundantCounter);
        System.out.println("Total runtime: " + ((double) (System.currentTimeMillis() - time) / 1000.0) + " seconds");
        // System.out.println(smth);
    }
}