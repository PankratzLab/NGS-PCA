package org.pankratzlab.ngspca;

import org.apache.commons.math3.linear.RealMatrix;
import org.ejml.data.DenseMatrix64F;
import org.genvisis.common.ArrayUtils;
import org.genvisis.stats.Maths;

/**
 * Common methods for operating on {@link RealMatrix}
 */
public class MatrixOperations {

	private MatrixOperations() {

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

	public enum SCALE_METHOD {
		/**
		 * compute fold-change (by column) , and then center the matrix so each row has
		 * median of 0;
		 */
		FC_MEDIAN,
		/**
		 * scale the columns of this matrix to a mean of 0 and SD of 1
		 */
		CENTER_SCALE_COLUMN,
		/**
		 * compute fold-change (by column) , and then center the matrix so each row has
		 * median of 0;
		 */
		CENTER_SCALE_COLUMN_SCALE_MARKER;
	}

	/**
	 * @param method
	 *            {@link SCALE_METHOD}
	 * @param m
	 *            {@link RealMatrix} to scale
	 */
	public static void scaleByMethod(SCALE_METHOD method, RealMatrix m) {
		switch (method) {
		case CENTER_SCALE_COLUMN:
			MatrixOperations.scaleAndCenterColumns(m);
			break;
		case CENTER_SCALE_COLUMN_SCALE_MARKER:
			MatrixOperations.scaleAndCenterColumns(m);
			MatrixOperations.centerRowsToMedian(m);
			break;
		case FC_MEDIAN:
			MatrixOperations.foldChangeAndCenter(m);
			break;
		default:
			throw new IllegalArgumentException("Invalid method " + method);
		}
	}

	/**
	 * compute fold-change (by column) , and then center the matrix so each row has
	 * median of 0;
	 * 
	 * @param m
	 *            an {@link RealMatrix} that has been FC-ed by column and centered
	 *            by row
	 */
	private static void foldChangeAndCenter(RealMatrix m) {
		double[] medians = new double[m.getColumnDimension()]; // In genvisis, samples are typically columns

		// convert columns to log2 fold-change from median
		for (int column = 0; column < m.getColumnDimension(); column++) {
			double[] tmp = new double[m.getRowDimension()];// In genvisis, markers are typically rows
			for (int row = 0; row < m.getRowDimension(); row++) {
				tmp[row] += m.getEntry(row, column);
			}
			medians[column] = ArrayUtils.median(tmp, true);
		}
		for (int row = 0; row < m.getRowDimension(); row++) {
			for (int column = 0; column < m.getColumnDimension(); column++) {
				double entry = m.getEntry(row, column);
				if (entry > 0) {
					double standard = Maths.log2(entry / medians[column]);
					m.setEntry(row, column, standard);
				} else {
					m.setEntry(row, column, 0);
				}
			}
		}
		// center rows to median of 0

		centerRowsToMedian(m);
	}

	/**
	 * @param m
	 *            Center the rows of this {@link RealMatrix} to a median of 0
	 */
	private static void centerRowsToMedian(RealMatrix m) {
		for (int row = 0; row < m.getRowDimension(); row++) {
			double[] tmp = m.getRow(row);
			double median = ArrayUtils.median(tmp, true);
			for (int column = 0; column < m.getColumnDimension(); column++) {
				m.setEntry(row, column, tmp[column] - median);
			}
		}
	}

	/**
	 * @param m
	 *            scale the columns of this matrix to a mean of 0 and SD of 1
	 */
	private static void scaleAndCenterColumns(RealMatrix m) {
		double[] sds = new double[m.getColumnDimension()];
		double[] mean = new double[m.getColumnDimension()];
		for (int column = 0; column < m.getColumnDimension(); column++) {
			double[] tmp = new double[m.getRowDimension()];
			for (int row = 0; row < m.getRowDimension(); row++) {
				tmp[row] += m.getEntry(row, column);
			}
			mean[column] = ArrayUtils.mean(tmp, true);
			sds[column] = ArrayUtils.stdev(tmp, true);
		}
		for (int row = 0; row < m.getRowDimension(); row++) {
			for (int column = 0; column < m.getColumnDimension(); column++) {
				double standard = m.getEntry(row, column) - mean[column];
				standard /= sds[column];
				m.setEntry(row, column, standard);
			}
		}
	}

}
