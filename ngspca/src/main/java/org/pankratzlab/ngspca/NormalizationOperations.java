package org.pankratzlab.ngspca;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.ranking.NaNStrategy;

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
  static void foldChangeAndCenterRows(RealMatrix dm) {
    // compute fold change
    computeFoldChangeByColumn(dm);
    // center rows to median of 0
    centerRowsToMedian(dm);
  }

  /**
   * Set the values of the matrix to the log 2 fold change (computed by column)
   * 
   * @param dm the {@link RealMatrix} that will be converted
   */
  private static void computeFoldChangeByColumn(RealMatrix dm) {
    double[] medians = new double[dm.getColumnDimension()];

    // convert columns to log2 fold-change from median
    for (int column = 0; column < dm.getColumnDimension(); column++) {
      double[] tmp = new double[dm.getRowDimension()];
      for (int row = 0; row < dm.getRowDimension(); row++) {
        tmp[row] += dm.getEntry(row, column);
      }
      medians[column] = median(tmp);
    }
    for (int row = 0; row < dm.getRowDimension(); row++) {
      for (int column = 0; column < dm.getColumnDimension(); column++) {
        double entry = dm.getEntry(row, column);
        if (entry > 0) {
          double standard = log2(entry / medians[column]);
          dm.setEntry(row, column, standard);
        } else {
          dm.setEntry(row, column, 0);
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
   * @param dm Center the rows of this {@link RealMatrix} to a median of 0
   */
  private static void centerRowsToMedian(RealMatrix dm) {
    for (int row = 0; row < dm.getRowDimension(); row++) {
      double[] tmp = new double[dm.getColumnDimension()];
      for (int col = 0; col < dm.getColumnDimension(); col++) {
        tmp[col] = dm.getEntry(row, col);
      }
      double median = median(tmp);
      for (int col = 0; col < dm.getColumnDimension(); col++) {
        dm.setEntry(row, col, tmp[col] - median);
      }
    }
  }

}
