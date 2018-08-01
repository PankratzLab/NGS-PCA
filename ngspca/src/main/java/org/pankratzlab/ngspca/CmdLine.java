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

/**
 * Utility class to process cmdline arguments
 *
 */
class CmdLine {
	static final String INPUT_DIR_ARG = "inputDir";
	static final String OUTPUT_DIR_ARG = "outputDir";
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
				.desc("print application usage and exit").hasArg(false).required(false).build();
		final Option inputOption = Option.builder("i").hasArg(true).longOpt(INPUT_DIR_ARG)
				.desc("Directory containing mosdepth results files (*.regions.bed.gz)").required().build();
		final Option outputOption = Option.builder("o").hasArg(true).required().longOpt(OUTPUT_DIR_ARG).hasArg()
				.desc("PCA results and auxillary files will be placed here").required().build();

		final Options options = new Options();
		options.addOption(inputOption);
		options.addOption(outputOption);
		options.addOption(help);

		return options;
	}

	/**
	 * Parse cmd line options
	 * 
	 * @param log
	 *            {@link Logger}
	 * @param options
	 *            {@link Options} to select from
	 * @param commandLineArguments
	 *            arguments provided to the cmd line
	 * @return parsed {@link CommandLine}
	 */
	static CommandLine generateCommandLine(Logger log, final Options options, final String[] commandLineArguments) {
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
	 * @param log
	 *            {@link Logger}
	 * @param options
	 *            {@link Options}
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