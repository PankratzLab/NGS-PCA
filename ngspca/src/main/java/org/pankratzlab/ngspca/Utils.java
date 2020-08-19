package org.pankratzlab.ngspca;

import java.util.logging.Level;
import java.util.logging.Logger;

public class Utils {

  static double[] convertToDoubleArray(String[] line, Logger log) {
    double[] values = new double[line.length];
    try {
      for (int i = 0; i < line.length; i++) {
        values[i] = Double.parseDouble(line[i]);
      }
    } catch (NumberFormatException nfe) {
      log.log(Level.SEVERE, "an exception was thrown", nfe);
      throw new IllegalArgumentException("Invalid (non-numeric) coverage value in matrix file");
    }
    return values;
  }
}
