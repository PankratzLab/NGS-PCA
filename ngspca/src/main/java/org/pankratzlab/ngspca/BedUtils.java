package org.pankratzlab.ngspca;

import java.io.Closeable;

import java.io.IOException;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;

public class BedUtils {

	private BedUtils() {

	}

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

	static BEDFileReader getReader(String file, boolean requireIndex) {
		return new BEDFileReader(file, requireIndex);
	}

	static class BEDFileReader implements Closeable, Iterable<BEDFeature> {

		private final FeatureReader<BEDFeature> reader;

		/**
		 * TODO : currently only handles .tbi indices
		 */
		private BEDFileReader(final String file, final boolean requireIndex) {
			reader = AbstractFeatureReader.getFeatureReader(file, new BEDCodec(), requireIndex);
		}

		// /** Queries for records within the region specified. */
		// public CloseableIterator<BEDFeature> query(final String chrom, final int
		// start, final int end) {
		// try {
		// return reader.query(chrom, start, end);
		// } catch (final IOException ioe) {
		// throw new TribbleException("Could not create an iterator from a feature
		// reader.", ioe);
		// }
		// }

		public void close() {
			try {
				reader.close();
			} catch (final IOException ioe) {
				throw new TribbleException("Could not close a bed context feature reader.", ioe);
			}
		}

		/** Returns an iterator over all records in this bed file */
		public CloseableIterator<BEDFeature> iterator() {
			try {
				return reader.iterator();
			} catch (final IOException ioe) {
				throw new TribbleException("Could not create an iterator from a feature reader.", ioe);
			}
		}

	}
}
