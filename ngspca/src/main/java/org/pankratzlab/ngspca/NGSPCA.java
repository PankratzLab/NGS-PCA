package org.pankratzlab.ngspca;

import java.io.File;
import java.util.HashSet;
import java.util.List;
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

  private static void run(String inputDir, String outputDir, boolean overwrite, Logger log) {
    new File(outputDir).mkdirs();

    String[] extensions = new String[] {MosdepthUtils.MOSDEPHT_BED_EXT};

    // get all files with mosdepth bed extension
    List<String> mosDepthResultFiles = FileOps.listFilesWithExtension(inputDir, extensions);
    mosDepthResultFiles = mosDepthResultFiles.subList(0, 30);

    if (mosDepthResultFiles.isEmpty()) {

      String err = "No input files provided";
      log.severe(err);
      throw new IllegalArgumentException(err);
    } else {
      log.info("Detected " + mosDepthResultFiles.size() + " mosdepth input files in " + inputDir
               + " with extension \"" + MosdepthUtils.MOSDEPHT_BED_EXT + "\"");
    }
    // parse sample names from files
    List<String> samples = mosDepthResultFiles.stream()
                                              .map(f -> FileOps.stripDirectoryAndExtension(f,
                                                                                           MosdepthUtils.MOSDEPHT_BED_EXT))
                                              .collect(Collectors.toList());

    // load ucsc regions to use

    List<String> regions = MosdepthUtils.getRegionsToUse(mosDepthResultFiles.get(0),
                                                         REGION_STRATEGY.AUTOSOMAL, log);

    String tmpDm = outputDir + "tmp.mat.ser.gz";

    String pcs = outputDir + "svd.pcs.txt";
    String loadings = outputDir + "svd.loadings.txt";
    String singularValues = outputDir + "svd.singularvalues.txt";

    DenseMatrix64F dm;
    if (!FileOps.fileExists(tmpDm)) {
      dm = MosdepthUtils.processFiles(mosDepthResultFiles, new HashSet<String>(regions), log);
    } else {
      dm = null;
    }

    SVD svd = new SVD(samples.toArray(new String[samples.size()]),
                      regions.toArray(new String[regions.size()]));
    svd.computeSVD(dm, 20, log);
    svd.dumpPCsToText(pcs, log);
    svd.computeAndDumpLoadings(loadings, dm, log);
    svd.dumpSingularValuesToText(singularValues, log);
  }

  public static void main(String[] args) {
    Logger log = Logger.getLogger(NGSPCA.class.getName());
    CommandLine cmd = CmdLine.generateCommandLine(log, CmdLine.generateOptions(), args);
    if (cmd == null || cmd.hasOption(CmdLine.HELP)) {
      CmdLine.printHelp(log, CmdLine.generateOptions());
      System.exit(1);
    }

    String inputDir = cmd.getOptionValue(CmdLine.INPUT_DIR_ARG);
    String outputDir = cmd.getOptionValue(CmdLine.OUTPUT_DIR_ARG);

    try {
      run(inputDir, outputDir, cmd.hasOption(CmdLine.OVERWRITE_ARG), log);
    } catch (Exception e) {
      log.log(Level.SEVERE, "an exception was thrown", e);
      log.severe("An exception occured while running\nFeel free to open an issue at https://github.com/PankratzLab/NGS-PCA after reviewing the help message below");
      CmdLine.printHelp(log, CmdLine.generateOptions());

    }
  }

  //
  // SerializedFiles.writeSerial(tm, svdFile, true);
  // }
  // log.reportTimeInfo("Loading " + svdFile);
  // SimpleNGSPCA tm = (SimpleNGSPCA) SerializedFiles.readSerial(svdFile, log,
  // false);
  // tm.dumpPCsToText(pcFile);
  // }
}
