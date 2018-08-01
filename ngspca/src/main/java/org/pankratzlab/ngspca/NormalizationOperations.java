package org.pankratzlab.ngspca;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.ejml.data.DenseMatrix64F;

/**
 * class to perform input matrix normalization prior to PCA
 */
public class NormalizationOperations {

	private NormalizationOperations() {

	}

	/**
	 * @param m
	 *            Convert this apache {@link RealMatrix} to an EJML style
	 *            {@link DenseMatrix64F}
	 * @return {@link DenseMatrix64F}
	 */
	public static DenseMatrix64F toDenseMatrix64F(RealMatrix m) {
		DenseMatrix64F dm = new DenseMatrix64F(1, 1);
		dm.reshape(m.getRowDimension(), m.getColumnDimension());
		for (int row = 0; row < m.getRowDimension(); row++) {
			for (int column = 0; column < m.getColumnDimension(); column++) {
				dm.add(row, column, m.getEntry(row, column));
			}
		}
		return dm;
	}

	/**
	 * compute fold-change (by column) , and then center the matrix so each row has
	 * median of 0;
	 * 
	 * @param m
	 *            an {@link RealMatrix} that has been FC-ed by column and centered
	 *            by row
	 */
	static void foldChangeAndCenterRows(RealMatrix m) {
		// compute fold change
		computeFoldChangeByColumn(m);
		// center rows to median of 0
		centerRowsToMedian(m);
	}

	/**
	 * Set the values of the matrix to the log 2 fold change (computed by column)
	 * 
	 * @param m
	 *            the {@link RealMatrix} that will be converted
	 */
	private static void computeFoldChangeByColumn(RealMatrix m) {
		double[] medians = new double[m.getColumnDimension()];

		// convert columns to log2 fold-change from median
		for (int column = 0; column < m.getColumnDimension(); column++) {
			double[] tmp = new double[m.getRowDimension()];
			for (int row = 0; row < m.getRowDimension(); row++) {
				tmp[row] += m.getEntry(row, column);
			}
			medians[column] = median(tmp);
		}
		for (int row = 0; row < m.getRowDimension(); row++) {
			for (int column = 0; column < m.getColumnDimension(); column++) {
				double entry = m.getEntry(row, column);
				if (entry > 0) {
					double standard = log2(entry / medians[column]);
					m.setEntry(row, column, standard);
				} else {
					m.setEntry(row, column, 0);
				}
			}
		}
	}

	private static double median(double[] tmp) {
		return new Median().withNaNStrategy(NaNStrategy.REMOVED).evaluate(tmp);
	}

	private static double log2(double num) {
		return Math.log(num) / Math.log(2);
	}

	/**
	 * @param m
	 *            Center the rows of this {@link RealMatrix} to a median of 0
	 */
	private static void centerRowsToMedian(RealMatrix m) {
		for (int row = 0; row < m.getRowDimension(); row++) {
			double[] tmp = m.getRow(row);
			double median = median(tmp);
			for (int column = 0; column < m.getColumnDimension(); column++) {
				m.setEntry(row, column, tmp[column] - median);
			}
		}
	}

}
