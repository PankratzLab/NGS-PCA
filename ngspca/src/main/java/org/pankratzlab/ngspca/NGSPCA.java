package org.pankratzlab.ngspca;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.pankratzlab.ngspca.BedUtils.BEDOverlapDetector;
import org.pankratzlab.ngspca.MosdepthUtils.REGION_STRATEGY;

/**
 * A simplified version of BamImport that uses MosDepth output or a custom input matrix to generate
 * PCS
 */
public class NGSPCA {

  private static void runInputMatrix(String inputMatrixFile, String outputDir, int numPcs,
                                     int niters, int numOversamples, int sampleAt, int randomSeed,
                                     boolean overwrite, boolean normMatrix,
                                     Logger log) throws InterruptedException, ExecutionException,
                                                 IOException {
    new File(outputDir).mkdirs();

    String delim = "\t";
    boolean gz = inputMatrixFile.endsWith(".gz");
    if (gz) {
      log.info("Assuming " + inputMatrixFile + " is gzipped");
    }

    log.info("Determining number of samples in " + inputMatrixFile);

    List<String> samples = FileOps.getFileHeader(inputMatrixFile, gz, delim, log);
    samples.remove(0);
    log.info("Found a total of " + samples.size() + " samples in " + inputMatrixFile);

    log.info("Determining number of regions in " + inputMatrixFile);
    List<String> regions = FileOps.getColumn(inputMatrixFile, gz, delim, 0, log);

    regions.remove(0);
    log.info("Found a total of " + regions.size() + " regions in " + inputMatrixFile);

    BlockRealMatrix dm;

    String tmpNormDm = outputDir + "tmp.mat.ser.gz";
    if (!FileOps.fileExists(tmpNormDm) || overwrite) {
      log.info("Populating matrix from " + inputMatrixFile);
      log.info("Initializing matrix to " + samples.size() + " columns and " + regions.size()
               + " rows");
      dm = new BlockRealMatrix(regions.size(), samples.size());
      int[] rowIndex = {0};
      //add data to matrix, skipping header and first column of file
      Stream<String> stream = gz ? FileOps.gzLines(Paths.get(inputMatrixFile), log)
                                 : Files.lines(Paths.get(inputMatrixFile));
      stream.skip(1).map(l -> l.split(delim))
            .forEach(a -> dm.setRow(rowIndex[0]++,
                                    Utils.convertToDoubleArray(Arrays.copyOfRange(a, 1, a.length),
                                                               0, log)));
      if (normMatrix) {
        log.info("Normalizing input matrix");
        NormalizationOperations.foldChangeAndCenterRows(dm, log);
      }
      FileOps.writeSerial(dm, tmpNormDm, log);
    } else {
      log.info("Loading existing serialized file " + tmpNormDm);
      dm = (BlockRealMatrix) FileOps.readSerial(tmpNormDm, log);
    }
    computeSVD(outputDir, numPcs, niters, numOversamples, randomSeed, log, samples, regions, dm);

  }

  /**
   * @param input directory or file listing full paths containing MosDepth results, with extension
   *          {@link MosdepthUtils#MOSDEPHT_BED_EXT}
   * @param outputDir where results will be written
   * @param bedExclude if not null, regions overlapping this bed file will not be included
   * @param regionStrategy how to select markers for PCA
   * @param numPcs number of PCs to retain in the output file
   * @param sampleAt sample the mosdepth bins, once per this number
   * @param randomSeed random seed for sampling matrix
   * @param overwrite overwrite any existing output
   * @param threads number of threads for loading bed files
   * @param log
   * @throws InterruptedException
   * @throws ExecutionException
   * @throws IOException
   */
  private static void runMosdepth(String input, String outputDir, String bedExclude,
                                  REGION_STRATEGY regionStrategy, int numPcs, int niters,
                                  int numOversamples, int sampleAt, int randomSeed,
                                  boolean overwrite, int threads,
                                  Logger log) throws InterruptedException, ExecutionException,
                                              IOException {
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
    // Store the raw input matrix
    String tmpRawDm = outputDir + "tmp.raw.ser.gz";
    // Store the temporary input matrix
    String tmpNormDm = outputDir + "tmp.mat.ser.gz";

    // populate input matrix and normalize
    BlockRealMatrix dm;
    if (!FileOps.fileExists(tmpNormDm) || overwrite) {
      dm = MosdepthUtils.processFiles(mosDepthResultFiles, new HashSet<String>(regions), tmpRawDm,
                                      threads, log);
      FileOps.writeSerial(dm, tmpNormDm, log);
    } else {
      System.out.print("Loading");
      System.err.print("Loading");
      log.info("Loading existing serialized file " + tmpNormDm);
      dm = (BlockRealMatrix) FileOps.readSerial(tmpNormDm, log);
    }
    //    String inputMatrix = outputDir + "svd.norm.input.txt";
    //    log.info("Writing to " + inputMatrix);
    //
    //    RandomizedSVD.dumpMatrix(inputMatrix, dm, "BIN", samples.toArray(new String[samples.size()]),
    //                             regions.toArray(new String[regions.size()]), false, log);

    computeSVD(outputDir, numPcs, niters, numOversamples, randomSeed, log, samples, regions, dm);
  }

  static void computeSVD(String outputDir, int numPcs, int niters, int numOversamples,
                         int randomSeed, Logger log, List<String> samples, List<String> regions,
                         BlockRealMatrix dm) {
    RandomizedSVD svd = new RandomizedSVD(samples, regions, log);

    log.info("Oversampling set to: " + numOversamples);
    log.info("Subspace iterations set to: " + niters);
    log.info("Random seed set to: " + randomSeed);
    svd.fit(dm, numPcs, niters, numOversamples, randomSeed);
    // perform SVD

    String pcs = outputDir + "svd.pcs.txt";
    String loadings = outputDir + "svd.loadings.txt";
    String singularValues = outputDir + "svd.singularvalues.txt";
    String binsUsed = outputDir + "svd.bins.txt";
    String samplesUsed = outputDir + "svd.samples.txt";

    log.info("Writing to " + pcs);
    svd.dumpPCsToText(pcs, log);
    log.info("Writing to " + loadings);
    svd.computeAndDumpLoadings(loadings, log);
    log.info("Writing to " + singularValues);
    svd.dumpSingularValuesToText(singularValues, log);
    log.info("Writing to " + binsUsed);
    FileOps.writeToText(svd.getRowNames(), binsUsed, log);
    log.info("Writing to " + samplesUsed);
    FileOps.writeToText(svd.getColumnNames(), samplesUsed, log);
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

      int niters = Integer.parseInt(cmd.getOptionValue(CmdLine.N_ITERS,
                                                       Integer.toString(RandomizedSVD.DEFAULT_NITERS)));
      int numOversamples = Integer.parseInt(cmd.getOptionValue(CmdLine.OVERSAMPLE,
                                                               Integer.toString(RandomizedSVD.DEFAULT_OVERSAMPLES)));

      int randomSeed = Integer.parseInt(cmd.getOptionValue(CmdLine.RANDOM_SEED,
                                                           Integer.toString(CmdLine.DEFAULT_RANDOM_SEED)));
      String bedExclude = cmd.getOptionValue(CmdLine.EXCLUDE_BED_FILE,
                                             CmdLine.DEFAULT_EXCLUDE_BED_FILE);
      if (cmd.hasOption(CmdLine.MATRIX_INPUT_ARG)) {
        runInputMatrix(input, outputDir, numPcs, niters, numOversamples, sampleAt, randomSeed,
                       cmd.hasOption(CmdLine.OVERWRITE_ARG),
                       cmd.hasOption(CmdLine.NORM_MATRIX_INPUT_ARG), log);
      } else {
        runMosdepth(input, outputDir, bedExclude, REGION_STRATEGY.AUTOSOMAL, numPcs, niters,
                    numOversamples, sampleAt, randomSeed, cmd.hasOption(CmdLine.OVERWRITE_ARG),
                    threads, log);
      }
    } catch (Exception e) {
      log.log(Level.SEVERE, "an exception was thrown", e);
      log.severe("An exception occured while running\nFeel free to open an issue at https://github.com/PankratzLab/NGS-PCA after reviewing the help message below");
      CmdLine.printHelp(log, CmdLine.generateOptions());
    }
  }
}
