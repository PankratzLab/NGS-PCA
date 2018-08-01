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

	/**
	 * @param inputFiles
	 *            mosdepth output bed files to be processed
	 * @param tmpFile
	 * @return
	 */
	static DenseMatrix64F processFiles(Logger log, List<String> inputFiles) {

		if (inputFiles.isEmpty()) {
			throw new IllegalArgumentException("No input files provided");
		}
		log.info("Starting input processing of " + inputFiles.size() + " files");

		Set<String> autosomalRegions = BedUtils.loadAutosomalUCSC(inputFiles.get(0));
		log.info("Initializing matrix to " + inputFiles.size() + " columns and " + autosomalRegions.size() + " rows");

		DenseMatrix64F dm = new DenseMatrix64F(1, 1);
		dm.reshape(autosomalRegions.size(), inputFiles.size());
		for (int col = 0; col < inputFiles.size(); col++) {
			if (col % 100 == 0) {
				log.info("Loading file " + Integer.toString(col + 1));
			}
			String inputFile = inputFiles.get(col);
			List<BEDFeature> features = BedUtils.loadSpecific(inputFile, autosomalRegions);
			for (int row = 0; row < features.size(); row++) {
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
