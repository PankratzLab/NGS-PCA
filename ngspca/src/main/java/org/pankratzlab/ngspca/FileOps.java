package org.pankratzlab.ngspca;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import java.util.logging.Level;
import java.util.logging.Logger;

public class FileOps {

  private FileOps() {

  }

  static String stripDirectoryAndExtension(String path, String extension) {

    return FilenameUtils.getName(path).replaceAll(extension, "");
  }

  static boolean fileExists(String path) {
    File f = new File(path);
    return f.exists() && !f.isDirectory();
  }

  static boolean dirExists(String path) {
    File f = new File(path);
    return f.isDirectory() && f.exists();
  }

  static boolean isDir(String path) {
    File f = new File(path);
    return f.isDirectory();
  }

  static List<String> readFile(String path) throws IOException {
    return FileUtils.readLines(new File(path), Charset.defaultCharset());
  }

  static List<String> listFilesWithExtension(String dir, String[] extensions) {
    return FileUtils.listFiles(new File(dir), extensions, true).stream().map(File::getAbsolutePath)
                    .collect(Collectors.toList());
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

  static void writeToText(List<String> list, String filename, Logger log) {

    try {
      FileWriter writer = new FileWriter(filename);
      for (String str : list) {
        writer.write(str + System.lineSeparator());

      }
      writer.close();
    } catch (IOException e1) {
      log.log(Level.SEVERE, "an exception was thrown while writing to " + filename, e1);
    }

  }

  static List<String> getColumn(String fileName, boolean gz, String delim, int columnNumber,
                                Logger log) {

    List<String> column = new ArrayList<>();

    try (Stream<String> stream = gz ? gzLines(Paths.get(fileName), log)
                                    : Files.lines(Paths.get(fileName))) {

      column = stream.map(s -> s.split(delim)[columnNumber])           // Split by ; and get the 2nd column
                     .collect(Collectors.toList());

    } catch (IOException e) {
      log.log(Level.SEVERE, "an exception was thrown while reading " + fileName, e);

    } catch (ArrayIndexOutOfBoundsException e) {
      log.log(Level.SEVERE, "invalid column index selected (" + columnNumber + ") in " + fileName,
              e);

    }
    return column;

  }

  static Stream<String> gzLines(Path path, Logger log) {
    InputStream fileIs = null;
    BufferedInputStream bufferedIs = null;
    GZIPInputStream gzipIs = null;
    try {
      fileIs = Files.newInputStream(path);
      bufferedIs = new BufferedInputStream(fileIs);
      gzipIs = new GZIPInputStream(bufferedIs);
    } catch (IOException e) {
      log.log(Level.SEVERE, "an exception was thrown while reading " + path.toString(), e);
      closeSafely(gzipIs, log);
      closeSafely(bufferedIs, log);
      closeSafely(fileIs, log);
    }
    BufferedReader reader = new BufferedReader(new InputStreamReader(gzipIs));
    return reader.lines().onClose(() -> closeSafely(reader, log));
  }

  private static void closeSafely(Closeable closeable, Logger log) {
    if (closeable != null) {
      try {
        closeable.close();
      } catch (IOException e) {
        log.log(Level.SEVERE, "an exception was thrown ", e);
      }
    }
  }

  static List<String> getFileHeader(String filename, boolean gz, String delim, Logger log) {
    List<String> header = new ArrayList<String>();
    try {
      Stream<String> stream = gz ? gzLines(Paths.get(filename), log)
                                 : Files.lines(Paths.get(filename));
      String[] tmpA = stream.findFirst().get().split(delim);
      for (String tmp : tmpA) {
        header.add(tmp);
      }
    } catch (FileNotFoundException e) {
      log.log(Level.SEVERE, "an exception was thrown while reading " + filename, e);
      return header;
    } catch (IOException e) {
      log.log(Level.SEVERE, "an exception was thrown while reading " + filename, e);

    }

    return header;
  }

}
