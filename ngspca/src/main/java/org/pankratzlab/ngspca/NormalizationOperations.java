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
   * compute fold-change (by column) , and then center the matrix so each row has median of 0;
   * 
   * @param m an {@link RealMatrix} that has been FC-ed by column and centered by row
   */
  static void foldChangeAndCenterRows(DenseMatrix64F dm) {
    // compute fold change
    computeFoldChangeByColumn(dm);
    // center rows to median of 0
    centerRowsToMedian(dm);
  }

  /**
   * Set the values of the matrix to the log 2 fold change (computed by column)
   * 
   * @param dm the {@link DenseMatrix64F} that will be converted
   */
  private static void computeFoldChangeByColumn(DenseMatrix64F dm) {
    double[] medians = new double[dm.getNumCols()];

    // convert columns to log2 fold-change from median
    for (int column = 0; column < dm.getNumCols(); column++) {
      double[] tmp = new double[dm.getNumRows()];
      for (int row = 0; row < dm.getNumRows(); row++) {
        tmp[row] += dm.get(row, column);
      }
      medians[column] = median(tmp);
    }
    for (int row = 0; row < dm.getNumRows(); row++) {
      for (int column = 0; column < dm.getNumCols(); column++) {
        double entry = dm.get(row, column);
        if (entry > 0) {
          double standard = log2(entry / medians[column]);
          dm.set(row, column, standard);
        } else {
          dm.set(row, column, 0);
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
   * @param dm Center the rows of this {@link DenseMatrix64F} to a median of 0
   */
  private static void centerRowsToMedian(DenseMatrix64F dm) {
    for (int row = 0; row < dm.getNumRows(); row++) {
      double[] tmp = new double[dm.getNumCols()];
      for (int col = 0; col < dm.getNumCols(); col++) {
        tmp[col] = dm.get(row, col);
      }
      double median = median(tmp);
      for (int col = 0; col < dm.getNumCols(); col++) {
        dm.set(row, col, tmp[col] - median);
      }
    }
  }

}
