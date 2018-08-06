package org.pankratzlab.ngspca;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.List;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import java.util.logging.Level;
import java.util.logging.Logger;

public class FileOps {

  static String stripDirectoryAndExtension(String path, String extension) {

    return FilenameUtils.getName(path).replaceAll(extension, "");
  }

  static boolean fileExists(String path) {
    File f = new File(path);
    return f.exists() && !f.isDirectory();
  }

  static List<String> listFilesWithExtension(String inputDir, String[] extensions) {
    return FileUtils.listFiles(new File(inputDir), extensions, true).stream()
                    .map(File::getAbsolutePath).collect(Collectors.toList());
  }

  static boolean writeSerial(Object o, String filename, Logger log) {

    try (ObjectOutputStream oos = new ObjectOutputStream(new GZIPOutputStream(new FileOutputStream(filename)))) {
      oos.writeObject(o);
      oos.flush();
      return true;
    } catch (Exception e) {
      log.log(Level.SEVERE, "an exception was thrown", e);
      return false;
    }
  }

  static Object readSerial(String filename, Logger log) {

    try (InputStream in = new BufferedInputStream(new GZIPInputStream(new FileInputStream(filename)))) {

      ObjectInputStream ois = new ObjectInputStream(in);
      return ois.readObject();
    } catch (Exception e) {
      log.log(Level.SEVERE, "an exception was thrown", e);
    }
    return null;

  }
}
