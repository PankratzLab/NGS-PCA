package org.pankratzlab.ngspca;

import java.util.List;
import java.util.Set;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingDeque;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.ejml.data.DenseMatrix64F;
import org.pankratzlab.ngspca.BedUtils.BEDOverlapDetector;
import org.pankratzlab.ngspca.BedUtils.BedRegionResult;
import htsjdk.tribble.bed.BEDFeature;

/**
 * Process mosdepth bed files for use in PCA
 */
class MosdepthUtils {

  static final String MOSDEPHT_BED_EXT = "regions.bed.gz";

  /**
   * TODO, region strategies for specific targets
   */
  enum REGION_STRATEGY {
    /**
     * Load all autosomal regions
     */
    AUTOSOMAL;
  }

  /**
   * @param mosDepthResultFile an example file to load regions from
   * @param rStrategy {@link REGION_STRATEGY} to use
   * @param excluder autosomal regions that return true for
   * @param log
   * @return {@link List} of regions
   */
  static List<String> getRegionsToUse(String mosDepthResultFile, REGION_STRATEGY rStrategy,
                                      BEDOverlapDetector excluder, Logger log) {

    log.info("Selecting regions using " + rStrategy + " region strategy");

    if (rStrategy == REGION_STRATEGY.AUTOSOMAL) {
      return BedUtils.loadAutosomalUCSC(mosDepthResultFile, excluder);
    } else {
      String err = "Invalid region strategy type " + rStrategy;
      log.severe(err);
      throw new IllegalArgumentException(err);
    }
  }

  /**
   * @param mosDepthResultFiles mosdepth output bed files to be processed
   * @param ucscRegions {@link Set} of regions to process
   * @param threads number of threads to use when loading
   * @param log
   * @return
   * @throws InterruptedException
   * @throws ExecutionException
   */
  static BlockRealMatrix processFiles(List<String> mosDepthResultFiles, Set<String> ucscRegions,
                                      String tmpRawFile, int threads,
                                      Logger log) throws InterruptedException, ExecutionException {
    if (mosDepthResultFiles.isEmpty()) {
      String err = "No input files provided";
      log.severe(err);
      throw new IllegalArgumentException(err);
    }
    return loadAndNormalizeData(mosDepthResultFiles, ucscRegions, tmpRawFile, threads, log);
  }

  /**
   * @param mosDepthResultFiles mosdepth output bed files to be processed
   * @param ucscRegions only these regions will be used
   * @param threads number of threads to use when loading
   * @param log
   * @return normalized {@link DenseMatrix64F} holding all input files
   * @throws InterruptedException
   * @throws ExecutionException
   */

  private static BlockRealMatrix loadAndNormalizeData(List<String> mosDepthResultFiles,
                                                      Set<String> ucscRegions, String tmpRawFile,
                                                      int threads, Logger log) {

    log.info("Initializing matrix to " + mosDepthResultFiles.size() + " columns and "
             + ucscRegions.size() + " rows");
    BlockRealMatrix dm = new BlockRealMatrix(ucscRegions.size(), mosDepthResultFiles.size());

    log.info("Starting input processing of " + mosDepthResultFiles.size() + " files");
    int col = 0;
    //    https://dzone.com/articles/the-evolution-of-producer-consumer-problem-in-java
    BlockingQueue<Future<BedRegionResult>> blockingQueue = new LinkedBlockingDeque<>(threads);
    ExecutorService executor = Executors.newFixedThreadPool(threads);

    Runnable producerTask = () -> {
      try {
        for (String mosDepthResultFile : mosDepthResultFiles) {
          blockingQueue.put(executor.submit(() -> BedUtils.loadSpecificRegions(mosDepthResultFile,
                                                                               ucscRegions)));
        }
      } catch (InterruptedException e) {
        log.severe(e.getMessage());
      }
    };

    executor.submit(producerTask);

    for (String mosDepthResultFile : mosDepthResultFiles) {
      try {
        BedRegionResult current = blockingQueue.take().get();
        String file = mosDepthResultFiles.get(col);
        if (!file.equals(current.file)) {
          throw new IllegalArgumentException("Invalid file returned, expecting " + file
                                             + " and got " + current.file);
        }
        if (current.features.size() != ucscRegions.size()) {
          throw new IllegalArgumentException("Invalid number of features from " + file
                                             + "\n expected" + ucscRegions.size() + " and got "
                                             + current.features.size());
        }
        setColumnData(dm, col, mosDepthResultFile, current.features, log);
        col++;
        if (col == 1 || col % 200 == 0) {
          log.info("Set data for file " + Integer.toString(col));
          log.info("Memory used: "
                   + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()));
        }
      } catch (InterruptedException e) {
        log.severe(e.getMessage());

      } catch (ExecutionException e) {
        log.severe(e.getMessage());

      }

    }
    executor.shutdown();
    log.info("Saving temporary raw matrix to " + tmpRawFile);
    FileOps.writeSerial(dm, tmpRawFile, log);

    log.info("Normalizing input matrix");
    NormalizationOperations.foldChangeAndCenterRows(dm, log);
    return dm;

  }

  private static void setColumnData(RealMatrix dm, int col, String inputFile,
                                    List<BEDFeature> features, Logger log) {

    for (int row = 0; row < features.size(); row++) {
      // mosdepth coverage parsed to "name" by htsjdk
      try {
        //        dm.add(v)
        //        dm.data
        dm.addToEntry(row, col, Double.parseDouble(features.get(row).getName()));
      } catch (NumberFormatException nfe) {
        log.log(Level.SEVERE, "an exception was thrown", nfe);
        throw new IllegalArgumentException("Invalid (non-numeric) coverage value in file "
                                           + inputFile + " in  row " + row);
      }

    }
  }

}
