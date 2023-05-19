package org.pankratzlab.ngspca;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.logging.Logger;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.pankratzlab.ngspca.RandomizedSVD.DISTRIBUTION;

/**
 * Utility class to process cmdline arguments
 */
class CmdLine {

  static final String INPUT_ARG = "input";
  static final String MATRIX_INPUT_ARG = "matrix";
  static final String NORM_MATRIX_INPUT_ARG = "normalizeMatrix";
  static final String OUTPUT_DIR_ARG = "outputDir";
  static final String OVERWRITE_ARG = "overwrite";
  static final String NUM_COMPONENTS_ARG = "numPC";
  static final String NUM_THREADS_ARG = "threads";
  static final String NUM_SAMPLE_ARG = "sampleEvery";
  static final String EXCLUDE_BED_FILE = "bedExclude";
  static final String N_ITERS = "iters";
  static final String OVERSAMPLE = "oversample";
  static final String RANDOM_SEED = "randomSeed";
  static final String DISTRIBUTION_ARG = "distribution";

  static final int DEFAULT_RANDOM_SEED = 42;
  static final int DEFAULT_PCS = 20;
  static final int DEFAULT_SAMPLE = 1;
  static final String DEFAULT_EXCLUDE_BED_FILE = null;
  static final DISTRIBUTION DEFAULT_DISTRIBUTION = DISTRIBUTION.UNIFORM;
  //  static final org.pankratzlab.ngspca.RandomizedSVD.DISTRIBUTION DEFAULT_DISTRIBUTION=DISTRIBUTION.
  //
  static final int DEFAULT_THREADS = 4;

  static final String HELP = "help";

  private CmdLine() {
    // utility
  }

  /**
   * Get available cmd line options
   * 
   * @return the {@link Options} available
   */
  static Options generateOptions() {
    final Option help = Option.builder("h").required().longOpt(HELP).hasArg()
                              .desc("Print application usage and exit").hasArg(false)
                              .required(false).build();
    final Option overwrite = Option.builder("r").required().longOpt(OVERWRITE_ARG).hasArg()
                                   .desc("Overwrite existing temporary files and recompute each step")
                                   .hasArg(false).required(false).build();
    final Option inputOption = Option.builder("i").hasArg(true).longOpt(INPUT_ARG)
                                     .desc("An existing directory containing mosdepth result files (*"
                                           + MosdepthUtils.MOSDEPHT_BED_EXT
                                           + " extension) OR a file listing paths to mosdepth result files, one result file per line OR a matrix of values to directly use (see the "
                                           + MATRIX_INPUT_ARG + " argument)")
                                     .required(true).build();

    final Option outputOption = Option.builder("o").hasArg(true).required().longOpt(OUTPUT_DIR_ARG)
                                      .hasArg()
                                      .desc("PCA results and auxillary files will be placed here")
                                      .required().build();
    final Option numComponents = Option.builder("n").hasArg(true).required()
                                       .longOpt(NUM_COMPONENTS_ARG).hasArg()
                                       .desc("The number of PCs to retain - the minimum of this and the number of markers/samples will be retained. Default is "
                                             + DEFAULT_PCS)
                                       .required(false).build();

    final Option numThreads = Option.builder("t").hasArg(true).required().longOpt(NUM_THREADS_ARG)
                                    .hasArg()
                                    .desc("Number of threads to utilize when loading data. Default is "
                                          + DEFAULT_THREADS)
                                    .required(false).build();

    final Option sampleEvery = Option.builder("s").hasArg(true).required().longOpt(NUM_SAMPLE_ARG)
                                     .hasArg()
                                     .desc("Sample mosdepth bins. For example, --" + NUM_SAMPLE_ARG
                                           + " 10 would select every tenth bin. Default is "
                                           + DEFAULT_SAMPLE + "(use every bin)")
                                     .required(false).build();

    final Option bedExcludes = Option.builder("b").hasArg(true).required().longOpt(EXCLUDE_BED_FILE)
                                     .hasArg()
                                     .desc("Optional: Provide a file to exclude specific regions from PCA input, prior to sampling with "
                                           + NUM_SAMPLE_ARG)
                                     .required(false).build();
    final Option niter = Option.builder(N_ITERS).hasArg(true).required().longOpt(N_ITERS).hasArg()
                               .desc("specifies the number of power (subspace) iterations to reduce the approximation error. The power scheme is recommended, if the singular values decay slowly. In practice, 2 or 3 iterations achieve good results, however, computing power iterations increases the computational costs "
                                     + RandomizedSVD.DEFAULT_NITERS)
                               .required(false).build();
    final Option oversamples = Option.builder(OVERSAMPLE).hasArg(true).required()
                                     .longOpt(OVERSAMPLE).hasArg()
                                     .desc("An oversampling parameter to improve the approximation of the randomized PCA. A value of at least 10 is recommended"
                                           + RandomizedSVD.DEFAULT_OVERSAMPLES)
                                     .required(false).build();

    final Option randomSeed = Option.builder(RANDOM_SEED).hasArg(true).required()
                                    .longOpt(RANDOM_SEED).hasArg()
                                    .desc("Random seed for generating sampling matrix for randomized PCA (probably not worth changing the default)"
                                          + DEFAULT_RANDOM_SEED)
                                    .required(false).build();

    final Option distribution = Option.builder(DISTRIBUTION_ARG).hasArg(true).required()
                                      .longOpt(DISTRIBUTION_ARG).hasArg()
                                      .desc("The distribution that will be used to seed the initial matrix. Options are "
                                            + RandomizedSVD.DISTRIBUTION.UNIFORM.toString() + " or "
                                            + RandomizedSVD.DISTRIBUTION.GAUSSIAN.toString() + "."
                                            + " Default is " + DEFAULT_DISTRIBUTION.toString())
                                      .required(false).build();
    final Option matrix = Option.builder(MATRIX_INPUT_ARG).hasArg(false).longOpt(MATRIX_INPUT_ARG)
                                .desc("The input provided by " + INPUT_ARG
                                      + " is a matrix (i.e. SVD will be performed directly on the matrix, without normalization, to generate PCS")
                                .required(false).build();
    final Option normMatrix = Option.builder(NORM_MATRIX_INPUT_ARG).hasArg(false)
                                    .longOpt(NORM_MATRIX_INPUT_ARG)
                                    .desc("The input provided by " + INPUT_ARG
                                          + " is a matrix and should be normalized (log2 by sample/column, centered by row)")
                                    .required(false).build();
    final Options options = new Options();
    options.addOption(help);

    options.addOption(inputOption);
    options.addOption(matrix);
    options.addOption(normMatrix);
    options.addOption(outputOption);
    options.addOption(numComponents);
    options.addOption(numThreads);
    options.addOption(sampleEvery);
    options.addOption(bedExcludes);
    options.addOption(niter);
    options.addOption(oversamples);
    options.addOption(randomSeed);
    options.addOption(distribution);
    options.addOption(overwrite);

    return options;
  }

  /**
   * Parse cmd line options
   * 
   * @param log {@link Logger}
   * @param options {@link Options} to select from
   * @param commandLineArguments arguments provided to the cmd line
   * @return parsed {@link CommandLine}
   */
  static CommandLine generateCommandLine(Logger log, final Options options,
                                         final String[] commandLineArguments) {
    final CommandLineParser cmdLineParser = new DefaultParser();
    CommandLine commandLine = null;
    try {
      commandLine = cmdLineParser.parse(options, commandLineArguments);
    } catch (ParseException parseException) {
      log.severe(parseException.getMessage());
    }
    return commandLine;
  }

  /**
   * Print the help usage to the log, based on {@link Options} available
   * 
   * @param log {@link Logger}
   * @param options {@link Options}
   */
  static void printHelp(Logger log, final Options options) {
    StringWriter out = new StringWriter();
    PrintWriter pw = new PrintWriter(out);

    HelpFormatter formatter = new HelpFormatter();
    formatter.printHelp(pw, 80, "ngspca", "USAGE:", options, 0, 0, "", true);
    pw.flush();
    pw.close();
    log.info(out.toString());
  }
}
