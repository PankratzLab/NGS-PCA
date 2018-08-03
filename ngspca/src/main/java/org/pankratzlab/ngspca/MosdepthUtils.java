package org.pankratzlab.ngspca;

import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

import org.ejml.data.DenseMatrix64F;

import htsjdk.tribble.bed.BEDFeature;

/**
 * Process mosdepth bed files for use in PCA
 *
 */
class MosdepthUtils {

	static final String MOSDEPHT_BED_EXT = "regions.bed.gz";

	/**
	 * TODO, region strategies for specific targets
	 *
	 */
	enum REGION_STRATEGY {
		/**
		 * Load all autosomal regions
		 */
		AUTOSOMAL;
	}

	/**
	 * @param mosDepthResultFile
	 *            an example file to load regions from
	 * @param rStrategy
	 *            {@link REGION_STRATEGY} to use
	 * @param log
	 * @return {@link List} of regions
	 */
	static List<String> getRegionsToUse(String mosDepthResultFile, REGION_STRATEGY rStrategy, Logger log) {

		switch (rStrategy) {
		case AUTOSOMAL:
			return BedUtils.loadAutosomalUCSC(mosDepthResultFile);
		default:
			String err = "Invalid region strategy type " + rStrategy;
			log.severe(err);
			throw new IllegalArgumentException(err);

		}
	}

	/**
	 * @param mosDepthResultFiles
	 *            mosdepth output bed files to be processed
	 * @param ucscRegions
	 *            {@link Set} of regions to process
	 * @param log
	 * @return
	 */
	static DenseMatrix64F processFiles(List<String> mosDepthResultFiles, Set<String> ucscRegions, Logger log) {
		if (mosDepthResultFiles.isEmpty()) {
			String err = "No input files provided";
			log.severe(err);
			throw new IllegalArgumentException(err);
		}
		return loadData(mosDepthResultFiles, ucscRegions, log);
	}

	/**
	 * @param mosDepthResultFiles
	 *            mosdepth output bed files to be processed
	 * @param ucscRegions
	 *            only these regions will be used
	 * @param log
	 * @return normalized {@link DenseMatrix64F} holding all input files
	 */

	private static DenseMatrix64F loadData(List<String> mosDepthResultFiles, Set<String> ucscRegions, Logger log) {

		// TODO use map to verify region indices
		log.info("Initializing matrix to " + mosDepthResultFiles.size() + " columns and " + ucscRegions.size()
				+ " rows");

		DenseMatrix64F dm = new DenseMatrix64F(1, 1);
		dm.reshape(ucscRegions.size(), mosDepthResultFiles.size());

		log.info("Starting input processing of " + mosDepthResultFiles.size() + " files");

		for (int col = 0; col < mosDepthResultFiles.size(); col++) {
			if (col % 100 == 0) {
				log.info("Loading file " + Integer.toString(col + 1));
			}
			String inputFile = mosDepthResultFiles.get(col);
			List<BEDFeature> features = BedUtils.loadSpecificRegions(inputFile, ucscRegions);
			for (int row = 0; row < features.size(); row++) {
				// mosdepth coverage parsed to "name" by htsjdk
				String tmp = features.get(row).getName();
				try {
					dm.set(row, col, Double.parseDouble(tmp));
				} catch (NumberFormatException nfe) {
					throw new IllegalArgumentException("Invalid value in file " + inputFile + ", row " + row);
				}

			}
		}
		log.info("Normalizing input matrix");
		NormalizationOperations.foldChangeAndCenterRows(dm);
		return dm;
	}

}
