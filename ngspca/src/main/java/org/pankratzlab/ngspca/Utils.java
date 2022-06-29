package org.pankratzlab.ngspca;

import java.util.logging.Level;
import java.util.logging.Logger;

public class Utils {

  static double[] convertToDoubleArray(String[] line, double defaultValue, Logger log) {
    double[] values = new double[line.length];

    for (int i = 0; i < line.length; i++) {
      try {
        values[i] = Double.parseDouble(line[i]);
      } catch (NumberFormatException nfe) {
        log.log(Level.SEVERE,
                "an exception was thrown on row " + i + " - setting to " + defaultValue, nfe);
        //        throw new IllegalArgumentException("Invalid (non-numeric) coverage value in matrix file");
        values[i] = defaultValue;
      }
    }

    return values;
  }
}
