package org.pankratzlab.ngspca;

import java.io.File;
import java.io.Serializable;
import java.util.HashSet;
import java.util.List;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.ejml.data.DenseMatrix64F;
import org.pankratzlab.ngspca.MosdepthUtils.REGION_STRATEGY;

/**
 * A simplified version of BamImport that uses MosDepth output to generate PCS
 */
public class SimpleNGSPCA implements Serializable {

	private static String parseSampleFromFilename(String path) {

		return FilenameUtils.getName(path).replaceAll(MosdepthUtils.MOSDEPHT_BED_EXT, "");
	}

	private static boolean fileExists(String path) {
		File f = new File(path);
		return f.exists() && !f.isDirectory();
	}

	private static void run(String inputDir, String outputDir, boolean overwrite, Logger log) {
		new File(outputDir).mkdirs();

		String[] extensions = new String[] { MosdepthUtils.MOSDEPHT_BED_EXT };

		// get all files with mosdepth bed extension
		List<String> mosDepthResultFiles = FileUtils.listFiles(new File(inputDir), extensions, true).stream()
				.map(File::getAbsolutePath).collect(Collectors.toList());
		if (mosDepthResultFiles.isEmpty()) {

			String err = "No input files provided";
			log.severe(err);
			throw new IllegalArgumentException(err);
		}
		// parse sample names from files

		List<String> samples = mosDepthResultFiles.stream().map(SimpleNGSPCA::parseSampleFromFilename)
				.collect(Collectors.toList());

		// load ucsc regions to use

		List<String> regions = MosdepthUtils.getRegionsToUse(mosDepthResultFiles.get(0), REGION_STRATEGY.AUTOSOMAL,
				log);

		String tmpDm = outputDir + "tmp.mat.ser.gz";
		DenseMatrix64F dm;
		if (!fileExists(tmpDm)) {
			dm = MosdepthUtils.processFiles(mosDepthResultFiles, new HashSet<String>(regions), log);
		} else {
			dm = null;
		}

		SVD svd = new SVD(samples.toArray(new String[samples.size()]), regions.toArray(new String[regions.size()]));
		svd.computeSVD(dm, log);
	}

	public static void main(String[] args) {
		Logger log = Logger.getLogger(SimpleNGSPCA.class.getName());
		CommandLine cmd = CmdLine.generateCommandLine(log, CmdLine.generateOptions(), args);
		if (cmd == null || cmd.hasOption(CmdLine.HELP)) {
			CmdLine.printHelp(log, CmdLine.generateOptions());
			System.exit(1);
		}

		String inputDir = cmd.getOptionValue(CmdLine.INPUT_DIR_ARG);
		String outputDir = cmd.getOptionValue(CmdLine.OUTPUT_DIR_ARG);
		run(inputDir, outputDir, cmd.hasOption(CmdLine.OVERWRITE_ARG), log);
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
