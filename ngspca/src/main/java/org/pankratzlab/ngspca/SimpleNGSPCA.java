package org.pankratzlab.ngspca;

import java.io.File;
import java.io.Serializable;
import java.util.List;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.io.FileUtils;

import org.pankratzlab.ngspca.MosdepthUtils.REGION_STRATEGY;

/**
 * A simplified version of BamImport that uses MosDepth output to generate PCS
 */
public class SimpleNGSPCA implements Serializable {

	private void dumpPCsToText(String file) {
		//
		// PrintWriter writer = Files.getAppropriateWriter(file);
		// StringJoiner joiner = new StringJoiner("\t");
		// joiner.add("SAMPLE");
		// for (int i = 0; i < v.getColumnDimension(); i++) {
		// joiner.add("PC" + (i + 1));
		// }
		// writer.println(joiner.toString());
		//
		// for (int i = 0; i < v.getRowDimension(); i++) {
		// StringJoiner sample = new StringJoiner("\t");
		// sample.add(colNames[i]);
		// for (int j = 0; j < v.getColumnDimension(); j++) {
		// sample.add(Double.toString(v.getEntry(j, i)));
		// }
		// writer.println(sample.toString());
		//
		// }
		// writer.close();
	}

	private static void run(String inputDir, String outputDir, Logger log) {
		new File(outputDir).mkdirs();

		String[] extensions = new String[] { MosdepthUtils.MOSDEPHT_BED_EXT };

		List<File> mosDepthResultFiles = (List<File>) FileUtils.listFiles(new File(inputDir), extensions, true);

		MosdepthUtils.processFiles(mosDepthResultFiles.stream().map(File::getAbsolutePath).collect(Collectors.toList()),
				REGION_STRATEGY.AUTOSOMAL, log);
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
		run(inputDir, outputDir, log);
	}

	// static void runGroup(CLI c, SCALE_METHOD method, String rootoutDir, Logger
	// log, String superPop,
	// Set<String> samps) {
	// String outDir = rootoutDir + superPop + "/" + method + "/";
	// new File(outDir).mkdirs();
	// log.reportTimeInfo("Operating in " + outDir);
	// String matFile = outDir + "mat.ser.gz";
	// String svdFile = outDir + "mat.ser.svd.gz";
	// String pcFile = outDir + "pcs.gz";
	//
	// if (!Files.exists(matFile)) {
	// List<String> files = new ArrayList<>();
	//
	// List<String> tmpFiles =
	// Arrays.asList(Files.listFullPaths(c.get(CLI.ARG_INDIR), "bed.gz"));
	// for (String tmpFile : tmpFiles) {
	// for (String sample : samps) {
	// if (tmpFile.contains(sample)) {
	// files.add(tmpFile);
	// break;
	// }
	// }
	// }
	//
	// log.reportTimeInfo("Found " + files.size() + " files to PCA for " +
	// superPop);
	// String[] availableRows = HashVec.loadFileToStringArray(files.get(0), false,
	// new int[] { 0, 1, 2 }, false);
	// List<String> rows = new ArrayList<>();
	// List<Integer> indices = new ArrayList<>();
	// for (int i = 0; i < availableRows.length; i++) {
	// String[] split = availableRows[i].trim().split("\t");
	//
	// Segment seg = new Segment(split[0], split[1], split[2]);
	// // or other filter logic
	// if (seg.getChr() > 0 && seg.getChr() < 23) {
	// rows.add(seg.getUCSClocation());
	// indices.add(i);
	// }
	// }
	// log.reportTimeInfo("Preparing matrix with " + rows.size() + " rows and " +
	// files.size() + "columns");
	//
	// RealMatrix m = MatrixUtils.createRealMatrix(rows.size(), files.size());
	// String[] colNames = new String[files.size()];
	// String[] rowNames = ArrayUtils.toStringArray(rows);
	// for (int i = 0; i < colNames.length; i++) {
	// colNames[i] = ext.rootOf(files.get(i));
	// log.reportTime("Parsing sample " + (i + 1) + ", " + colNames[i]);
	// String[] data = HashVec.loadFileToStringArray(files.get(i), false, new int[]
	// { 0, 1, 2, 3 }, false);
	// for (int j = 0; j < indices.size(); j++) {
	// String[] current = data[indices.get(j)].split("\t");
	// Segment seg = new Segment(current[0], current[1], current[2]);
	// if (!seg.getUCSClocation().equals(rows.get(j))) {
	// throw new IllegalStateException(
	// "Should be " + rows.get(j) + "and found " + seg.getUCSClocation());
	// }
	//
	// m.addToEntry(j, i, Double.parseDouble(current[3]));
	// }
	// }
	// // new SingularValueDecomposition(
	//
	// log.reportTime("Saving progress");
	//
	// SerializedFiles.writeSerial(new SimpleNGSPCA(m, colNames, rowNames), matFile,
	// true);
	// }
	// if (!Files.exists(svdFile)) {
	// log.reportTime("Loading data");
	// SimpleNGSPCA tm = (SimpleNGSPCA) SerializedFiles.readSerial(matFile, log,
	// false);
	// log.reportTime("Scaling data");
	// MatrixOperations.scaleByMethod(method, tm.m);
	//
	// log.reportTimeInfo("Rows=" + tm.m.getRowDimension());
	// log.reportTimeInfo("Cols=" + tm.m.getColumnDimension());
	// log.reportTime("Computing SVD");
	// tm.computeSVD(log);
	//
	// log.reportTime("Finished Computing SVD");
	//
	// log.reportTime("Saving progress");
	//
	// SerializedFiles.writeSerial(tm, svdFile, true);
	// }
	// log.reportTimeInfo("Loading " + svdFile);
	// SimpleNGSPCA tm = (SimpleNGSPCA) SerializedFiles.readSerial(svdFile, log,
	// false);
	// tm.dumpPCsToText(pcFile);
	// }
}
