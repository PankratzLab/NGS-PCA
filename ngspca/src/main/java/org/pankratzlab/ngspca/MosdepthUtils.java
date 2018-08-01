package org.pankratzlab.ngspca;

import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

import org.apache.commons.lang3.StringUtils;
import org.ejml.data.DenseMatrix64F;

import htsjdk.tribble.bed.BEDFeature;

/**
 * Process mosdepth bed files for use in PCA
 *
 */
class MosdepthUtils {

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

	// = BedUtils.loadAutosomalUCSC(mosDepthResultFiles.get(0));

	/**
	 * @param mosDepthResultFiles
	 *            mosdepth output bed files to be processed
	 * @param ucscRegions
	 *            only these regions will be used
	 * @param log
	 * @return normalized {@link DenseMatrix64F} holding all input files
	 */

	private static DenseMatrix64F processFiles(List<String> mosDepthResultFiles, Set<String> ucscRegions, Logger log) {

		if (mosDepthResultFiles.isEmpty()) {
			throw new IllegalArgumentException("No input files provided");
		}

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
				if (!StringUtils.isNumeric(tmp)) {
					throw new IllegalArgumentException("Invalid value in file " + inputFile + ", row " + row);
				}
				dm.set(row, col, Double.parseDouble(tmp));

			}
		}
		log.info("Normalizing input matrix");
		NormalizationOperations.foldChangeAndCenterRows(dm);
		return dm;
	}

}
