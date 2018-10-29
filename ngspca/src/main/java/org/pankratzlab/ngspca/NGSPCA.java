package org.pankratzlab.ngspca;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import org.apache.commons.cli.CommandLine;
import org.ejml.data.DenseMatrix64F;
import org.pankratzlab.ngspca.BedUtils.BEDOverlapDetector;
import org.pankratzlab.ngspca.MosdepthUtils.REGION_STRATEGY;

/**
 * A simplified version of BamImport that uses MosDepth output to generate PCS
 */
public class NGSPCA {

  /**
   * @param input directory or file listing full paths containing MosDepth results, with extension
   *          {@link MosdepthUtils#MOSDEPHT_BED_EXT}
   * @param outputDir where results will be written
   * @param bedExclude if not null, regions overlapping this bed file will not be included
   * @param regionStrategy how to select markers for PCA
   * @param numPcs number of PCs to retain in the output file
   * @param sampleAt sample the mosdepth bins, once per this number
   * @param overwrite overwrite any existing output
   * @param threads number of threads for loading bed files
   * @param log
   * @throws InterruptedException
   * @throws ExecutionException
   * @throws IOException
   */
  private static void run(String input, String outputDir, String bedExclude,
                          REGION_STRATEGY regionStrategy, int numPcs, int sampleAt,
                          boolean overwrite, int threads,
                          Logger log) throws InterruptedException, ExecutionException, IOException {
    new File(outputDir).mkdirs();

    String[] extensions = new String[] {MosdepthUtils.MOSDEPHT_BED_EXT};

    // get all files with mosdepth bed extension
    List<String> mosDepthResultFiles;
    if (FileOps.isDir(input) && FileOps.dirExists(input)) {
      log.info("Detected " + input + " is a directory, searching for "
               + MosdepthUtils.MOSDEPHT_BED_EXT + " extensions");
      mosDepthResultFiles = FileOps.listFilesWithExtension(input, extensions);
    } else if (FileOps.fileExists(input)) {
      log.info("Detected " + input + " is a file, reading mosdepth result file paths");
      mosDepthResultFiles = FileOps.readFile(input);
    } else {
      String err = "Invalid or non-existent input argument " + input;
      log.severe(err);
      throw new IllegalArgumentException(err);
    }

    if (mosDepthResultFiles.isEmpty()) {

      String err = "No input files were found";
      log.severe(err);
      throw new IllegalArgumentException(err);
    } else {
      log.info("Detected " + mosDepthResultFiles.size() + " mosdepth input files in " + input);
    }
    // parse sample names from files
    List<String> samples = mosDepthResultFiles.stream()
                                              .map(f -> FileOps.stripDirectoryAndExtension(f,
                                                                                           MosdepthUtils.MOSDEPHT_BED_EXT))
                                              .collect(Collectors.toList());

    // load ucsc regions to use

    BEDOverlapDetector overlapDetector = new BEDOverlapDetector(bedExclude, log);
    List<String> regions = MosdepthUtils.getRegionsToUse(mosDepthResultFiles.get(0), regionStrategy,
                                                         overlapDetector, log);
    log.info(overlapDetector.getNumExcluded() + " regions removed during up-front filtering");
    if (sampleAt > 1) {
      log.info("Sampling the" + regions.size() + " mosdepth regions once every " + sampleAt
               + " bins");
      regions = IntStream.range(0, regions.size()).filter(n -> n % sampleAt == 0)
                         .mapToObj(regions::get).collect(Collectors.toList());
      log.info("Sampled " + regions.size() + " bins");

    }
    String tmpDm = outputDir + "tmp.mat.ser.gz";

    // populate input matrix and normalize
    DenseMatrix64F dm;
    if (!FileOps.fileExists(tmpDm) || overwrite) {
      dm = MosdepthUtils.processFiles(mosDepthResultFiles, new HashSet<String>(regions), threads,
                                      log);
      FileOps.writeSerial(dm, tmpDm, log);
    } else {
      log.info("Loading existing serialized file " + tmpDm);
      dm = (DenseMatrix64F) FileOps.readSerial(tmpDm, log);
    }

    // perform SVD
    SVD svd = new SVD(samples.toArray(new String[samples.size()]),
                      regions.toArray(new String[regions.size()]));
    svd.computeSVD(dm, numPcs, log);

    // output results
    String pcs = outputDir + "svd.pcs.txt";
    String loadings = outputDir + "svd.loadings.txt";
    String singularValues = outputDir + "svd.singularvalues.txt";

    log.info("Writing to " + pcs);
    svd.dumpPCsToText(pcs, log);
    log.info("Writing to " + loadings);
    svd.computeAndDumpLoadings(loadings, dm, log);
    log.info("Writing to " + singularValues);
    svd.dumpSingularValuesToText(singularValues, log);
  }

  public static void main(String[] args) {
    Logger log = Logger.getLogger(NGSPCA.class.getName());
    CommandLine cmd = CmdLine.generateCommandLine(log, CmdLine.generateOptions(), args);
    if (cmd == null || cmd.hasOption(CmdLine.HELP)) {
      CmdLine.printHelp(log, CmdLine.generateOptions());
      System.exit(1);
    }

    String input = cmd.getOptionValue(CmdLine.INPUT_ARG);
    String outputDir = cmd.getOptionValue(CmdLine.OUTPUT_DIR_ARG);

    try {
      int numPcs = Integer.parseInt(cmd.getOptionValue(CmdLine.NUM_COMPONENTS_ARG,
                                                       Integer.toString(CmdLine.DEFAULT_PCS)));

      int threads = Integer.parseInt(cmd.getOptionValue(CmdLine.NUM_THREADS_ARG,
                                                        Integer.toString(CmdLine.DEFAULT_THREADS)));

      int sampleAt = Integer.parseInt(cmd.getOptionValue(CmdLine.NUM_SAMPLE_ARG,
                                                         Integer.toString(CmdLine.DEFAULT_SAMPLE)));

      String bedExclude = cmd.getOptionValue(CmdLine.EXCLUDE_BED_FILE,
                                             CmdLine.DEFAULT_EXCLUDE_BED_FILE);
      run(input, outputDir, bedExclude, REGION_STRATEGY.AUTOSOMAL, numPcs, sampleAt,
          cmd.hasOption(CmdLine.OVERWRITE_ARG), threads, log);
    } catch (Exception e) {
      log.log(Level.SEVERE, "an exception was thrown", e);
      log.severe("An exception occured while running\nFeel free to open an issue at https://github.com/PankratzLab/NGS-PCA after reviewing the help message below");
      CmdLine.printHelp(log, CmdLine.generateOptions());

    }
  }
}
