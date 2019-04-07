package org.pankratzlab.ngspca;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.linear.MatrixUtils;
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
  static RealMatrix processFiles(List<String> mosDepthResultFiles, Set<String> ucscRegions,
                                 int threads,
                                 Logger log) throws InterruptedException, ExecutionException {
    if (mosDepthResultFiles.isEmpty()) {
      String err = "No input files provided";
      log.severe(err);
      throw new IllegalArgumentException(err);
    }
    return loadAndNormalizeData(mosDepthResultFiles, ucscRegions, threads, log);
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

  private static RealMatrix loadAndNormalizeDataOld(List<String> mosDepthResultFiles,
                                                    Set<String> ucscRegions, int threads,
                                                    Logger log) throws InterruptedException,
                                                                ExecutionException {

    // TODO use map to verify region indices
    log.info("Initializing matrix to " + mosDepthResultFiles.size() + " columns and "
             + ucscRegions.size() + " rows");

    //    RealMatrix dm = new (ucscRegions.size(), mosDepthResultFiles.size());
    RealMatrix dm = MatrixUtils.createRealMatrix(ucscRegions.size(), mosDepthResultFiles.size());
    //        new DenseMatrix64F(ucscRegions.size(), mosDepthResultFiles.size());

    log.info("Starting input processing of " + mosDepthResultFiles.size() + " files");

    ExecutorService executor = Executors.newFixedThreadPool(threads);
    List<Future<BedRegionResult>> futures = new ArrayList<>();
    for (String mosDepthResultFile : mosDepthResultFiles) {
      futures.add(executor.submit(() -> BedUtils.loadSpecificRegions(mosDepthResultFile,
                                                                     ucscRegions)));
    }

    int col = 0;
    for (Future<BedRegionResult> future : futures) {
      BedRegionResult current = future.get();
      //      future.
      String file = mosDepthResultFiles.get(col);
      if (!file.equals(current.file)) {
        throw new IllegalArgumentException("Invalid file returned, expecting " + file + " and got "
                                           + current.file);
      }
      setColumnData(dm, col, mosDepthResultFiles.get(col), current.features, log);
      col++;

      if (col == 1 || col % 100 == 0) {
        log.info("Set data for file " + Integer.toString(col));
        log.info("Memory used: "
                 + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()));
      }
    }
    log.info("Normalizing input matrix");
    NormalizationOperations.foldChangeAndCenterRows(dm);
    executor.shutdown();
    return dm;

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

  private static RealMatrix loadAndNormalizeData(List<String> mosDepthResultFiles,
                                                 Set<String> ucscRegions, int threads, Logger log) {

    log.info("Initializing matrix to " + mosDepthResultFiles.size() + " columns and "
             + ucscRegions.size() + " rows");
    RealMatrix dm = MatrixUtils.createRealMatrix(ucscRegions.size(), mosDepthResultFiles.size());

    log.info("Starting input processing of " + mosDepthResultFiles.size() + " files");
    int col = 0;
    for (String mosDepthResultFile : mosDepthResultFiles) {

      setColumnData(dm, col, mosDepthResultFile,
                    BedUtils.loadSpecificRegions(mosDepthResultFile, ucscRegions).features, log);
      col++;
      if (col == 1 || col % 200 == 0) {
        log.info("Set data for file " + Integer.toString(col));
        log.info("Memory used: "
                 + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()));
      }
    }
    log.info("Normalizing input matrix");
    NormalizationOperations.foldChangeAndCenterRows(dm);
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
