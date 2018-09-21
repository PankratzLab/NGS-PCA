package org.pankratzlab.ngspca;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ExecutionException;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import org.apache.commons.cli.CommandLine;
import org.ejml.data.DenseMatrix64F;
import org.pankratzlab.ngspca.MosdepthUtils.REGION_STRATEGY;

/**
 * A simplified version of BamImport that uses MosDepth output to generate PCS
 */
public class NGSPCA {

  /**
   * @param input directory or file listing full paths containing MosDepth results, with extension
   *          {@link MosdepthUtils#MOSDEPHT_BED_EXT}
   * @param outputDir where results will be written
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
  private static void run(String input, String outputDir, REGION_STRATEGY regionStrategy,
                          int numPcs, int sampleAt, boolean overwrite, int threads,
                          Logger log) throws InterruptedException, ExecutionException, IOException {
    new File(outputDir).mkdirs();

    String[] extensions = new String[] {MosdepthUtils.MOSDEPHT_BED_EXT};

    // get all files with mosdepth bed extension
    List<String> mosDepthResultFiles = new ArrayList<>();
    if (FileOps.isDir(input) && FileOps.dirExists(input)) {
      log.info("Detected " + input + " is a directory, searching for "
               + MosdepthUtils.MOSDEPHT_BED_EXT + " extensions");
      mosDepthResultFiles = FileOps.listFilesWithExtension(input, extensions);
    } else if (FileOps.fileExists(input)) {
      log.info("Detected " + input + " is a file, reading mosdepth result file paths");
      mosDepthResultFiles = FileOps.readFile(input);
    } else {
      String err = "Invalid or non-existent input argument ";
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

    List<String> regions = MosdepthUtils.getRegionsToUse(mosDepthResultFiles.get(0), regionStrategy,
                                                         log);
    if (sampleAt > 1) {
      log.info("Sampling the" + regions.size() + " mosdepth regions once every " + sampleAt
               + " bins");
      List<String> tmp = new ArrayList<>();

      for (int i = 0; i < regions.size(); i++) {
        if (i % sampleAt == 0) {
          tmp.add(regions.get(i));
        }
      }
      regions = tmp;
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
      run(input, outputDir, REGION_STRATEGY.AUTOSOMAL, numPcs, sampleAt,
          cmd.hasOption(CmdLine.OVERWRITE_ARG), threads, log);
    } catch (Exception e) {
      log.log(Level.SEVERE, "an exception was thrown", e);
      log.severe("An exception occured while running\nFeel free to open an issue at https://github.com/PankratzLab/NGS-PCA after reviewing the help message below");
      CmdLine.printHelp(log, CmdLine.generateOptions());

    }
  }
}
