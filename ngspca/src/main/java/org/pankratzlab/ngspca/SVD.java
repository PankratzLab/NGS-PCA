package org.pankratzlab.ngspca;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;
import java.util.StringJoiner;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.io.FileUtils;
import org.ejml.alg.dense.decomposition.svd.SvdImplicitQrDecompose_D64;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import org.ejml.ops.SingularOps;

/**
 * class to perform singular value decomposition on {@link DenseMatrix64F}
 */
class SVD implements Serializable {

  /**
   * 
   */
  private static final long serialVersionUID = 5517545544997813383L;
  /**
   * Column names of the original input data
   */
  private final List<String> originalColNames;
  /**
   * Row names of the original input data
   */
  private final List<String> originalRowNames;

  private SvdImplicitQrDecompose_D64 svd;

  /**
   * Number of components that are stored
   */
  private int numComponents;

  /**
   * @param colNames
   * @param rowNames
   */
  SVD(List<String> colNames, List<String> rowNames) {
    super();
    this.originalColNames = colNames;
    this.originalRowNames = rowNames;
  }

  /**
   * @param dm {@link DenseMatrix64F} that has been pre-normalized, if desired
   * @param numberOfComponentsToStore Number of components (PCs) that are stored
   * @param log
   */

  void computeSVD(DenseMatrix64F dm, int numberOfComponentsToStore, Logger log) {
    if (dm.getNumCols() != originalColNames.size()) {
      throw new IllegalArgumentException("Mismatched column lengths");
    }
    if (dm.getNumRows() != originalRowNames.size()) {
      throw new IllegalArgumentException("Mismatched row lengths");
    }
    log.info("Computing EJML PCs");
    this.svd = new SvdImplicitQrDecompose_D64(false, false, true, true);
    svd.decompose(dm);
    log.info("Finished Computing EJML PCs");

    //    svd.getU(null, false)

    log.info("Sorting matrices to descending order");

    SingularOps.descendingOrder(null, false, svd.getW(null), svd.getV(null, true), true);

    this.numComponents = Math.min(numberOfComponentsToStore,
                                  Math.min(svd.getW(null).numRows, svd.getV(null, true).numCols));
    if (numComponents < numberOfComponentsToStore) {
      log.info(numberOfComponentsToStore + " PCs requested, but only be able to compute "
               + numComponents);
    }
    log.info("Subsetting matrices to " + numComponents + " components");

    svd.getV(null, true).reshape(numComponents, dm.getNumCols(), true);

  }

  /**
   * @param numComponents number of singular values to return from the W matrix
   * @return
   */
  double[] getSingularValues() {
    double[] singularValues = new double[numComponents];
    for (int i = 0; i < numComponents; i++) {
      singularValues[i] = svd.getW(null).get(i, i);
    }
    return singularValues;
  }

  /**
   * @param file dump the PCs to this text file
   * @param log
   */
  void dumpPCsToText(String file, Logger log) {
    //
    DenseMatrix64F v = svd.getV(null, true);

    List<String> pcNames = getNumberedColumnHeader("PC", v.getNumRows());

    dumpMatrix(file, v, "SAMPLE", pcNames, originalColNames, true, log);
  }

  /**
   * @param file loadings will be computed and dumped to this file
   * @param dm the original {@link DenseMatrix64F}
   * @param log
   */
  void computeAndDumpLoadings(String file, DenseMatrix64F dm, Logger log) {
    DenseMatrix64F loadingData = computeLoadings(dm);
    List<String> loadingNames = getNumberedColumnHeader("Loading", loadingData.getNumCols());
    dumpMatrix(file, loadingData, "MARKER", loadingNames, originalRowNames, false, log);

  }

  /**
   * @param file singular values will be dumped to this file
   * @param log
   */
  void dumpSingularValuesToText(String file, Logger log) {
    StringJoiner joiner = new StringJoiner("\n");
    joiner.add("SINGULAR_VALUES\tPC");

    for (int component = 0; component < numComponents; component++) {
      joiner.add(component + "\t" + Double.toString(svd.getW(null).get(component, component)));
    }

    try {
      FileUtils.writeStringToFile(new File(file), joiner.toString(), Charset.defaultCharset(),
                                  false);
    } catch (IOException e) {
      log.log(Level.SEVERE, "unable to write to file " + file, e);
    }

  }

  public static List<String> getNumberedColumnHeader(String type, int num) {
    List<String> names = new ArrayList<>();
    for (int i = 0; i < num; i++) {
      names.add(type + (i + 1));
    }
    return names;
  }

  private static void dumpMatrix(String file, DenseMatrix64F m, String rowTitle,
                                 List<String> outputColumnNames, List<String> outputRowNames,
                                 boolean transposed, Logger log) {
    try (PrintWriter writer = new PrintWriter(new File(file))) {

      StringJoiner joiner = new StringJoiner("\t");
      joiner.add(rowTitle);
      for (String colName : outputColumnNames) {
        joiner.add(colName);
      }
      writer.println(joiner.toString());
      log.info(outputRowNames.size() + " output rows by " + outputColumnNames.size()
               + " output columns");

      for (int outputRow = 0; outputRow < outputRowNames.size(); outputRow++) {
        StringJoiner sample = new StringJoiner("\t");
        sample.add(outputRowNames.get(outputRow));
        for (int outputColumn = 0; outputColumn < outputColumnNames.size(); outputColumn++) {

          //          sample.add(Double.toString(m.get(j, i)));
          if (transposed) {
            sample.add(Double.toString(m.get(outputColumn, outputRow)));
          } else {
            sample.add(Double.toString(m.get(outputRow, outputColumn)));
          }

        }
        writer.println(sample.toString());

      }
    } catch (FileNotFoundException e) {

      log.log(Level.SEVERE, "unable to write to file " + file, e);

    }
  }

  private DenseMatrix64F computeLoadings(DenseMatrix64F dm) {
    //    Will have all markers, but not all "PCs" all the time
    DenseMatrix64F loadingData = new DenseMatrix64F(dm.numRows, numComponents);
    double[] singularValues = getSingularValues();
    for (int row = 0; row < dm.numRows; row++) {

      DenseMatrix64F rowData = new DenseMatrix64F(1, dm.numCols);
      CommonOps.extractRow(dm, row, rowData);

      for (int component = 0; component < numComponents; component++) {
        DenseMatrix64F componentData = new DenseMatrix64F(1, svd.getV(null, true).numCols);

        double loading = getLoading(singularValues[component], rowData.data,
                                    CommonOps.extractRow(svd.getV(null, true), component,
                                                         componentData).data);
        loadingData.add(row, component, loading);
      }
    }
    return loadingData;
  }

  private static double getLoading(double singularValue, double[] data, double[] basis) {
    double loading = 0;
    for (int i = 0; i < basis.length; i++) {
      loading += data[i] * basis[i];
    }
    return loading / singularValue;
  }

  /**
   * @param serFile write to this serialized file
   */
  void writeSerial(String serFile, Logger log) {
    FileOps.writeSerial(this, serFile, log);
  }

  /**
   * @param serFile
   * @return {@link SVD}
   */
  static SVD readSerial(String serFile, Logger log) {
    return (SVD) FileOps.readSerial(serFile, log);
  }

}
