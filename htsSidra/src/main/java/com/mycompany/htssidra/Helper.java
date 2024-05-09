package com.mycompany.htssidra;

import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.Reader;
import java.io.Writer;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import net.sourceforge.htmlunit.corejs.javascript.ast.DoLoop;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;
import org.apache.commons.lang3.ArrayUtils;

/**
 *
 * @author ealiyev
 */
public class Helper {

    public static void calculateSVTypeSvaba(String filename) throws FileNotFoundException, IOException {
        try {
            CSVReader reader = new CSVReader(new FileReader(filename), '\t');
            List<String[]> rows = reader.readAll();
            String[] array = new String[]{"##INFO=<ID=END,Number=1,Type=Integer,Description=\"END position\">"};
            rows.add(109, array);
            for (Iterator<String[]> iterator = rows.iterator(); iterator.hasNext();) {
                String[] row = iterator.next();
                if (!row[0].startsWith("#")) {
                    try {
                        String info = row[7];
                        String altValue = row[4];
                        String idValue = row[2];
                        String chrom = row[0];
                        String start = row[1];
                        String mateChrom = "";
                        mateChrom = altValue.replaceAll("[^0-9:XY]", "").split(":")[0];
                        if (info.contains("SECONDARY")
                                || idValue.contains(":2")
                                || chrom.startsWith("GL")
                                || chrom.contains("hs37d5")
                                || chrom.contains("MT")
                                || altValue.contains("hs37d5")
                                || altValue.contains("GL")
                                || !chrom.equals(mateChrom)) {
                            iterator.remove();
                            continue;
                        }
                        String svInfo = "";
                        if (altValue.indexOf("]") > 0) {
                            svInfo = "++";
                            row[7] = row[7].replace("SVTYPE=BND", "SVTYPE=INV");
                        } else if (altValue.indexOf("]") == 0) {
                            svInfo = "-+";
                            row[7] = row[7].replace("SVTYPE=BND", "SVTYPE=DUP");
                        } else if (altValue.indexOf("[") > 0) {
                            svInfo = "+-";
                            row[7] = row[7].replace("SVTYPE=BND", "SVTYPE=DEL");
                        } else if (altValue.indexOf("[") == 0) {
                            svInfo = "--";
                            row[7] = row[7].replace("SVTYPE=BND", "SVTYPE=INV");
                        }
                        if (row[7].contains("SPAN")) {
                            row[7] = row[7].replace("SPAN", "SVLEN");
                        }

                        String[] infoValues = row[7].split(";");
                        for (int i = 0; i < infoValues.length; i++) {
                            if (infoValues[i].split("=").length > 1 && infoValues[i].split("=")[0].equals("SVLEN")) {
                                Integer svLength = Integer.parseInt(infoValues[i].split("=")[1]);
                                Integer startPos = Integer.valueOf(start);
                                Integer endPos = Integer.sum(svLength, startPos);
                                row[7] = row[7] + ";END=" + String.valueOf(endPos) + ";";
                            }
                        }

                    } catch (Exception e) {
                        System.out.println("Error in format");
                        e.printStackTrace();
                    }

                }

            }

            String replace = filename.replace(".svaba.sv.vcf", ".converted.vcf");
            CSVWriter writer = new CSVWriter(new FileWriter(replace), '\t', CSVWriter.NO_QUOTE_CHARACTER);
            writer.writeAll(rows);
            Thread.sleep(2000);
            writer.close();
        } catch (InterruptedException ex) {
            Logger.getLogger(Helper.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    private static boolean criteriaFilter(HashMap<String, String> infoValuesMap) {

        String chromValue = infoValuesMap.get("CHROM");
        String altValue = infoValuesMap.get("ALT");
        if (chromValue.startsWith("GL") || chromValue.contains("hs37d5") || chromValue.contains("MT") || altValue.contains("hs37") || altValue.contains("GL")) {
            return false;
        }
        return true;
    }

    static void mergeFiles(String arg) throws IOException, InterruptedException {
        File[] vcfList = finder(arg, ".tsv");
        List<String[]> finalRows = new ArrayList<>();
        for (int i = 0; i < vcfList.length; i++) {
            List<String[]> tempRows = Helper.parsevcfINFO(vcfList[i]);
            if (i == 0) {
                finalRows.addAll(tempRows);
            } else {
                if (tempRows.size() > 0) {
                    tempRows.remove(0);
                    finalRows.addAll(tempRows);
                }
            }
        }
        saveMergedFile(arg, finalRows);

    }

    public static File[] finder(String dirName, String extension) {
        File dir = new File(dirName);
        return dir.listFiles(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String filename) {
                return filename.endsWith(extension);
            }
        });
    }

    public static void saveMergedFile(String path, List<String[]> rows) {
        FileWriter fileWriter = null;
        CSVPrinter csvFilePrinter = null;
        try {
            fileWriter = new FileWriter(path + "merged.tsv");

            CSVFormat csvFileFormat = CSVFormat.MYSQL;

            csvFilePrinter = new CSVPrinter(fileWriter, csvFileFormat);

            csvFilePrinter.printRecords(rows);

        } catch (IOException ex) {
            Logger.getLogger(Helper.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                fileWriter.flush();
                fileWriter.close();
                csvFilePrinter.close();
            } catch (IOException ex) {
                Logger.getLogger(Helper.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    public static List<String[]> readAnnotatedVCF(File file) {
        try {
            CSVReader reader = new CSVReader(new FileReader(file), '\t');
            List<String[]> rows = reader.readAll();
            return rows;
        } catch (FileNotFoundException ex) {
            Logger.getLogger(Helper.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(Helper.class.getName()).log(Level.SEVERE, null, ex);
        }
        return null;
    }

    public static List<String[]> parsevcfINFO(File file) throws IOException, InterruptedException {
        FileWriter fileWriter = null;
        CSVPrinter csvFilePrinter = null;

        List<String[]> rows = Helper.readAnnotatedVCF(file);
        int infoColumnPos = 0;
        for (int i = 0; i < rows.size(); i++) {
            String[] row = rows.get(i);
            if (!(row[0].contains("CHROM") || row[0].startsWith("#"))) {
                String[] split = row[infoColumnPos].split(";");

                HashMap<String, String> info = new HashMap<>();
                for (String split1 : split) {
                    if (split1.split("=").length > 1) {
                        info.put(split1.split("=")[0], split1.split("=")[1]);
                    }
                }
                Integer svLen = null;
                Integer svEnd = null;
                boolean isBig = false;
                if (!info.getOrDefault("SVLEN", "0").equals("0")) {
                    svLen = Math.abs(Integer.valueOf(info.getOrDefault("SVLEN", "0")));
                    svEnd = Math.abs(Integer.valueOf(info.getOrDefault("END", "0")));
                }

                List<String> list = new ArrayList<>(Arrays.asList(row));
                list.add(2, String.valueOf(svLen));
                list.add(2, String.valueOf(svEnd));
                list.add(0, "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr" + row[0] + "%3A" + row[1] + "-" + svEnd + "&hgsid=602254539_aQwveDR5zDGOQeFrQTPJjkhX0LqF");
                list.add(0, file.getName());

                row = list.toArray(new String[0]);

                rows.set(i, row);

            } else if (row[0].contains("CHROM")) {
                for (int j = 0; j < row.length; j++) {
                    if (row[j].equals("INFO")) {
                        infoColumnPos = j;
                    }
                }
                List<String> list = new ArrayList<>(Arrays.asList(row));

                list.add(2, "SVLEN");
                list.add(2, "END");
                list.add(0, "USCS_URL");
                list.add(0, "Family");

                row = list.toArray(new String[0]);

                rows.set(i, row);

            }

            fileWriter = new FileWriter("C:\\Users\\ealiyev\\Dropbox\\SV Project\\AnnotatedPMC\\PMC_HOMO_16_08_2017\\" + file.getName().replace("annot", "annot_trans"));

            CSVFormat csvFileFormat = CSVFormat.MYSQL;

            csvFilePrinter = new CSVPrinter(fileWriter, csvFileFormat);

        }
        System.out.println(file);
        System.out.println(rows.size());
        csvFilePrinter.printRecords(rows);

        fileWriter.flush();
        fileWriter.close();
        csvFilePrinter.close();
        return rows;
//       
    }

    static void popVCFtoPoppante(String filename) throws FileNotFoundException, IOException {
        try {
            CSVReader reader = new CSVReader(new FileReader(filename), '\t');
            List<String[]> rows = reader.readAll();
            List<String> sampleList = new ArrayList<>();
            List<String[]> idData = new ArrayList<>();

            List<String[]> bamFiles = new ArrayList<>();

            for (Iterator<String[]> iterator = rows.iterator(); iterator.hasNext();) {
                String[] row = iterator.next();
                if (row[0].startsWith("#CHROM")) {
                    try {
                        for (int i = 9; i < row.length; i++) {
                            sampleList.add(row[i]);
                        }
                        System.out.println(sampleList);
                    } catch (Exception e) {
                        System.out.println("Error in format");
                        e.printStackTrace();
                    }

                } else if (!row[0].contains("#")) {
                    idData.clear();
                    String currentId = row[2];

                    for (int i = 9; i < row.length; i++) {
                        String[] cur_format = row[i].split(":");
                        switch (cur_format[0]) {
                            case "1/1":
                                idData.add(new String[]{sampleList.get(i - 9).split("_")[0], sampleList.get(i - 9).split("_")[1], "2"});
                                break;
                            case "0/1":
                                idData.add(new String[]{sampleList.get(i - 9).split("_")[0], sampleList.get(i - 9).split("_")[1], "1"});
                                break;
                            default:
                                idData.add(new String[]{sampleList.get(i - 9).split("_")[0], sampleList.get(i - 9).split("_")[1], "0"});
                                break;
                        }
                    }

                    CSVWriter writer = new CSVWriter(new FileWriter(currentId + ".predictor"), '\t', CSVWriter.NO_QUOTE_CHARACTER);
                    writer.writeAll(idData);
                    Thread.sleep(1000);
                    writer.close();

                    PrintWriter mapwriter = new PrintWriter(currentId + ".map", "UTF-8");
                    mapwriter.println(currentId);
                    mapwriter.close();
                }

            }

        } catch (InterruptedException ex) {
            Logger.getLogger(Helper.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    static void renameSampleName(String arg) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    static void extractSampleBam(String arg) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    static void processAnnotatedFile(String filename) {
        try {
            CSVReader reader = new CSVReader(new FileReader(filename), '\t');
            List<String[]> rows = reader.readAll();
            ArrayList<String[]> processedRows = new ArrayList<>();
            CSVWriter writer = new CSVWriter(new FileWriter(filename.trim() + ".processed.tsv"), '\t', CSVWriter.NO_QUOTE_CHARACTER);

            HashMap<String, Integer> headerHelper = new HashMap<>();

            HashMap<String, String[]> rowbyId = new HashMap();

            for (Iterator<String[]> iterator = rows.iterator(); iterator.hasNext();) {
                String[] next = iterator.next();
                Boolean novel = false;
                Boolean genic = false;
                Boolean intronic = false;
                Boolean exonic = false;
                Boolean promoterSite = false;
                Boolean trascriptionSite = false;
                Boolean pLOF = false;
                Boolean DiseaseGenes = false;
                double highestPLI = 0.0;
                if (next[0].contains("SV")) {
                    for (int i = 0; i < next.length; i++) {
                        headerHelper.put(next[i], i);
                    }
                    String[] added = {"exon_merged", "pli_merged", "DDD_disease_merged", "omim_numbers", "GenePanels_merged", "phenotypes_merged", 
                        "highest_PLI", "Genic", "Exonic", "Intronic", "Novel", "Promoter", "pLOF", "DISEASE_GENES"};
                    String[] newRow = (String[]) ArrayUtils.addAll(next, added);
                    writer.writeNext(newRow);
                } else if (next[headerHelper.get("AnnotSV type")].equals("full")) {
                    List<String> location = new ArrayList<>();
                    List<Float> pli = new ArrayList<>();
                    List<String> phenotypes = new ArrayList<>();
                    List<String> dddPhenotypes = new ArrayList<>();
                    List<String> omim_numbers = new ArrayList<>();
                    List<String> genepanels = new ArrayList<>();
                    for (String[] row : rows) {
                        if (next[headerHelper.get("ID")].equals(row[headerHelper.get("ID")]) && row[headerHelper.get("AnnotSV type")].equals("split")) {
                            location.add(row[headerHelper.get("location")]);
                            dddPhenotypes.add(row[headerHelper.get("DDD_disease")]);
                            omim_numbers.add(row[headerHelper.get("Mim Number")]);
                            if (headerHelper.getOrDefault("GenePanel", 0) != 0) {
                                genepanels.add(row[headerHelper.get("GenePanel")]);
                            }
                            try {
                                pli.add(Float.valueOf(row[headerHelper.get("pLI_ExAC")]));
                            } catch (Exception e) {

                            }
                            phenotypes.add(row[headerHelper.get("Phenotypes")]);
                        }

                    }

                    DiseaseGenes = Pattern.compile("[0-9]").matcher(phenotypes.toString()).find();

                    if (pli.size() > 0) {
                        highestPLI = Collections.max(pli);
                    }
                    genic = !next[headerHelper.get("Gene name")].isEmpty();
                    promoterSite = !next[headerHelper.get("promoters")].isEmpty();
                    exonic = !next[headerHelper.get("Gene name")].isEmpty() && (location.toString().contains("tx") || location.toString().contains("ex"));
                    if (next[headerHelper.get("SV type")].equals("DEL") && exonic) {
                        pLOF = true;
                    }
                    for (int i = 0; i < location.size(); i++) {
                        String curLoc = location.get(i);
                        if ((next[headerHelper.get("SV type")].equals("DUP") || next[headerHelper.get("SV type")].equals("INV")) && curLoc.contains("txStart-txEnd")) {
                            pLOF = true;
                        }
                        if (curLoc.contains("intron")) {
                            String[] locs = curLoc.split("-");
                            if (!locs[0].replace("intron", "").equals(locs[1].replace("intron", ""))) {
                                exonic = true;
                                if (next[headerHelper.get("SV type")].equals("DEL")) {
                                    pLOF = true;
                                }
                            }
                        }
                    }

                    intronic = !exonic;

                    novel = next[headerHelper.get("DGV_GAIN_IDs")].isEmpty()
                            && next[headerHelper.get("DGV_LOSS_IDs")].isEmpty()
                            && next[headerHelper.get("DDD_SV")].isEmpty()
                            && next[headerHelper.get("1000g_event")].isEmpty()
                            && next[headerHelper.get("GD_ID")].isEmpty()
                            && next[headerHelper.get("IMH_ID")].isEmpty()
                            && next[headerHelper.get("dbVar_variant")].isEmpty();

                    String[] added = {
                        Arrays.toString(location.toArray()).replace("[", "").replace("]", ""),
                        Arrays.toString(pli.toArray()).replace("[", "").replace("]", ""),
                        Arrays.toString(dddPhenotypes.toArray()).replace("[", "").replace("]", ""),
                        Arrays.toString(omim_numbers.toArray()).replace("[", "").replace("]", ""),
                        genepanels.toString(),
                        Arrays.toString(phenotypes.toArray()).replace("[", "").replace("]", ""),
                        String.valueOf(highestPLI),
                        genic.toString(),
                        exonic.toString(),
                        intronic.toString(),
                        novel.toString(),
                        promoterSite.toString(),
                        pLOF.toString(),
                        DiseaseGenes.toString()};
                    String[] newRow = (String[]) ArrayUtils.addAll(next, added);
                    writer.writeNext(newRow);
                }

            }
            writer.close();

        } catch (IOException ex) {
            Logger.getLogger(Helper.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    static void processAnnotatedFile_3_8(String filename) {
        System.out.println("Processing annotated file");
        try {
            CSVReader reader = new CSVReader(new FileReader(filename), '\t');
            List<String[]> rows = reader.readAll();
            CSVWriter writer = new CSVWriter(new FileWriter(filename.trim() + ".processed.tsv"), '\t', CSVWriter.NO_QUOTE_CHARACTER);

            HashMap<String, Integer> headerHelper = new HashMap<>();

            HashMap<String, String[]> rowbyId = new HashMap();

            for (Iterator<String[]> iterator = rows.iterator(); iterator.hasNext();) {
                String[] next = iterator.next();
                Boolean novel = false;
                Boolean genic = false;
                Boolean intronic = false;
                Boolean exonic = false;
                Boolean promoterSite = false;
                Boolean pLOF = false;
                Boolean DiseaseGenes = false;
                double highestPLI = 0.0;
                if (next[0].contains("SV")) {
                    for (int i = 0; i < next.length; i++) {
                        headerHelper.put(next[i], i);
                    }

                    String[] added = {"tx_merged","frameshift_merged", 
                        "Exon_count_merged", "Overlapped_CDS_length_merged", 
                        "Location_merged", "Location2_merged", 
                        "Dist_nearest_SS_merged", "Nearest_SS_type_merged", 
                        "ACMG_merged", "HI_merged", "TS_merged","DDD_HI_PERCENT",
                        "GenCC_disease_merged","GenCC_moi_merged",
                        "GenCC_classification_merged","GenCC_pmid_merged","OMIM_ID_merged",
                        "OMIM_phenotype_merged","OMIM_inheritance_merged","OMIM_morbid_merged",
                        "OMIM_morbid_candidate_merged","LOEUF_bin_merged",
                        "highest_PLI", "Genic", "Exonic", "Intronic", "Novel",
                        "pLOF", "DISEASE_GENES"};
                    String[] newRow = (String[]) ArrayUtils.addAll(next, added);
                    writer.writeNext(newRow);
                } else if (next[headerHelper.get("Annotation_mode")].equals("full")) {
                    List<String> tx_merged = new ArrayList<>();
                    List<String> frameshift_merged = new ArrayList<>();
                    List<String> Exon_count_merged = new ArrayList<>();
                    List<String> Overlapped_CDS_length_merged = new ArrayList<>();
                    List<String> Location_merged = new ArrayList<>();
                    List<String> Location2_merged = new ArrayList<>();
                    List<String> Dist_nearest_SS_merged = new ArrayList<>();
                    List<String> Nearest_SS_type_merged = new ArrayList<>();
                    List<String> ACMG_merged = new ArrayList<>();
                    List<Float> HI_merged = new ArrayList<>();
                    List<Float> TS_merged = new ArrayList<>();
                    List<Float> DDD_HI_PERCENT_merged = new ArrayList<>();
                    List<String> GenCC_disease_merged = new ArrayList<>();
                    List<String> GenCC_moi_merged = new ArrayList<>();
                    List<String> GenCC_classification_merged = new ArrayList<>();
                    List<String> GenCC_pmid_merged = new ArrayList<>();
                    List<String> OMIM_ID_merged = new ArrayList<>();
                    List<String> OMIM_phenotype_merged = new ArrayList<>();
                    List<String> OMIM_inheritance_merged = new ArrayList<>();
                    List<String> OMIM_morbid_merged = new ArrayList<>();
                    List<String> OMIM_morbid_candidate_merged = new ArrayList<>();
                    List<Float> LOEUF_bin_merged = new ArrayList<>();                           
                    
                    List<Float> pli = new ArrayList<>();

                    for (String[] row : rows) {
                        if (next[headerHelper.get("ID")].equals(row[headerHelper.get("ID")]) && row[headerHelper.get("Annotation_mode")].equals("split")) {
                           tx_merged.add(row[headerHelper.get("Tx")]);
                           frameshift_merged.add(row[headerHelper.get("Frameshift")]);
                           Exon_count_merged.add(row[headerHelper.get("Exon_count")]);
                           Overlapped_CDS_length_merged.add(row[headerHelper.get("Overlapped_CDS_length")]);
                           Location_merged.add(row[headerHelper.get("Location")]);
                           Location2_merged.add(row[headerHelper.get("Location2")]);
                           Dist_nearest_SS_merged.add(row[headerHelper.get("Dist_nearest_SS")]);
                           Nearest_SS_type_merged.add(row[headerHelper.get("Nearest_SS_type")]);
                           ACMG_merged.add(row[headerHelper.get("ACMG")]);
                           GenCC_disease_merged.add(row[headerHelper.get("GenCC_disease")]);
                           GenCC_moi_merged.add(row[headerHelper.get("GenCC_moi")]);
                           GenCC_classification_merged.add(row[headerHelper.get("GenCC_classification")]);
                           GenCC_pmid_merged.add(row[headerHelper.get("GenCC_pmid")]);
                           OMIM_ID_merged.add(row[headerHelper.get("OMIM_ID")]);
                           OMIM_phenotype_merged.add(row[headerHelper.get("OMIM_phenotype")]);
                           OMIM_inheritance_merged.add(row[headerHelper.get("OMIM_inheritance")]);
                           OMIM_morbid_merged.add(row[headerHelper.get("OMIM_morbid")]);
                           OMIM_morbid_candidate_merged.add(row[headerHelper.get("OMIM_morbid_candidate")]);

    
                            try {
                                pli.add(Float.valueOf(row[headerHelper.get("pLI_ExAC")]));
                                HI_merged.add(Float.valueOf(row[headerHelper.get("HI")]));
                                TS_merged.add(Float.valueOf(row[headerHelper.get("TS")]));
                                DDD_HI_PERCENT_merged.add(Float.valueOf(row[headerHelper.get("DDD_HI_percent")]));
                                LOEUF_bin_merged.add(Float.valueOf(row[headerHelper.get("LOEUF_bin")]));
                            } catch (Exception e) {
                            }
                        }

                    }

                    DiseaseGenes = Pattern.compile("[0-9]").matcher(OMIM_phenotype_merged.toString()).find();

                    if (pli.size() > 0) {
                        highestPLI = Collections.max(pli);
                    }
                    genic = !next[headerHelper.get("Gene_name")].isEmpty();
                    exonic = !next[headerHelper.get("Gene_name")].isEmpty() && (Location_merged.toString().contains("tx") || Location_merged.toString().contains("ex"));
                    if (next[headerHelper.get("SV_type")].equals("DEL") && exonic) {
                        pLOF = true;
                    }
                    for (int i = 0; i < Location_merged.size(); i++) {
                        String curLoc = Location_merged.get(i);
                        if ((next[headerHelper.get("SV_type")].equals("DUP") || next[headerHelper.get("SV_type")].equals("INV")) && curLoc.contains("txStart-txEnd")) {
                            pLOF = true;
                        }
                        if (curLoc.contains("intron")) {
                            String[] locs = curLoc.split("-");
                            if (!locs[0].replace("intron", "").equals(locs[1].replace("intron", ""))) {
                                exonic = true;
                                if (next[headerHelper.get("SV_type")].equals("DEL")) {
                                    pLOF = true;
                                }
                            }
                        }
                    }

                    intronic = !exonic;

                    novel = false;
                    
                    String[] added = {
                        Arrays.toString(tx_merged.toArray()).replace("[", "").replace("]", ""),
                        Arrays.toString(frameshift_merged.toArray()).replace("[", "").replace("]", ""),
                        Arrays.toString(Exon_count_merged.toArray()).replace("[", "").replace("]", ""),
                        Arrays.toString(Overlapped_CDS_length_merged.toArray()).replace("[", "").replace("]", ""),
                        Arrays.toString(Location_merged.toArray()).replace("[", "").replace("]", ""),
                        Arrays.toString(Location2_merged.toArray()).replace("[", "").replace("]", ""),
                        Arrays.toString(Dist_nearest_SS_merged.toArray()).replace("[", "").replace("]", ""),
                        Arrays.toString(Nearest_SS_type_merged.toArray()).replace("[", "").replace("]", ""),
                        Arrays.toString(ACMG_merged.toArray()).replace("[", "").replace("]", ""),
                        Arrays.toString(HI_merged.toArray()).replace("[", "").replace("]", ""),
                        Arrays.toString(TS_merged.toArray()).replace("[", "").replace("]", ""),
                        Arrays.toString(DDD_HI_PERCENT_merged.toArray()).replace("[", "").replace("]", ""),
                        Arrays.toString(GenCC_disease_merged.toArray()).replace("[", "").replace("]", ""),
                        Arrays.toString(GenCC_moi_merged.toArray()).replace("[", "").replace("]", ""),
                        Arrays.toString(GenCC_classification_merged.toArray()).replace("[", "").replace("]", ""),
                        Arrays.toString(GenCC_pmid_merged.toArray()).replace("[", "").replace("]", ""),
                        Arrays.toString(OMIM_ID_merged.toArray()).replace("[", "").replace("]", ""),
                        Arrays.toString(OMIM_phenotype_merged.toArray()).replace("[", "").replace("]", ""),
                        Arrays.toString(OMIM_inheritance_merged.toArray()).replace("[", "").replace("]", ""),
                        Arrays.toString(OMIM_morbid_merged.toArray()).replace("[", "").replace("]", ""),
                        Arrays.toString(OMIM_morbid_candidate_merged.toArray()).replace("[", "").replace("]", ""),
                        Arrays.toString(LOEUF_bin_merged.toArray()).replace("[", "").replace("]", ""),
                        String.valueOf(highestPLI),
                        genic.toString(),
                        exonic.toString(),
                        intronic.toString(),
                        novel.toString(),
                        pLOF.toString(),
                        DiseaseGenes.toString()};
                        
                        String[] newRow = (String[]) ArrayUtils.addAll(next, added);
                        writer.writeNext(newRow);
                }

            }
            writer.close();

        } catch (IOException ex) {
            Logger.getLogger(Helper.class.getName()).log(Level.SEVERE, null, ex);
        }
    }


    static void processAnnotatedFileSplit(String filename) {
        
        try {
            CSVReader reader = new CSVReader(new FileReader(filename), '\t');
            List<String[]> rows = reader.readAll();
            ArrayList<String[]> processedRows = new ArrayList<>();
            CSVWriter writer = new CSVWriter(new FileWriter(filename + ".processed.tsv"), '\t', CSVWriter.NO_QUOTE_CHARACTER);

            HashMap<String, Integer> headerHelper = new HashMap<>();

            Boolean novel = false;
            Boolean genic = false;
            Boolean intronic = false;
            Boolean exonic = false;
            Boolean promoterSite = false;

            HashMap<String, String[]> rowbyId = new HashMap();

            for (Iterator<String[]> iterator = rows.iterator(); iterator.hasNext();) {
                String[] next = iterator.next();

                if (next[0].contains("SV")) {
                    for (int i = 0; i < next.length; i++) {
                        headerHelper.put(next[i], i);
                    }
                } else if (next[headerHelper.get("AnnotSV type")].equals("full")) {

                } else if (next[headerHelper.get("AnnotSV type")].equals("split")) {

                    genic = !next[headerHelper.get("Gene name")].isEmpty();
                    promoterSite = !next[headerHelper.get("promoters")].isEmpty();
                    exonic = next[headerHelper.get("location")].contains("ex") || next[headerHelper.get("location")].contains("tx");
                    String curLoc = next[headerHelper.get("location")];
                    if (curLoc.contains("intron")) {
                        String[] locs = curLoc.split("-");
                        if (!locs[0].replace("intron", "").equals(locs[1].replace("intron", ""))) {
                            exonic = true;
                        }
                    }

                    intronic = !exonic;

                    novel = next[headerHelper.get("DGV_GAIN_IDs")].isEmpty()
                            && next[headerHelper.get("DGV_LOSS_IDs")].isEmpty()
                            && next[headerHelper.get("DDD_SV")].isEmpty()
                            && next[headerHelper.get("1000g_event")].isEmpty() && next[headerHelper.get("dbVar_variant")].isEmpty();

                    String[] added = {genic.toString(), exonic.toString(), intronic.toString(), novel.toString(), promoterSite.toString()};
                    String[] newRow = (String[]) ArrayUtils.addAll(next, added);
                    writer.writeNext(newRow);
                }

            }
            writer.close();

        } catch (IOException ex) {
            Logger.getLogger(Helper.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    static void vcfStats(String filename, String annotation, String filter) {
        try {
            CSVReader reader = new CSVReader(new FileReader(filename), '\t');
            CSVReader annotations = new CSVReader(new FileReader(annotation), '\t');

            List<String[]> rows = reader.readAll();
            List<String[]> annotationrows = annotations.readAll();
            HashMap<String, Integer> headerHelper = new HashMap<>();
            HashMap<String, String[]> rowbyID = new HashMap<>();

            List<String> sampleList = new ArrayList<>();

            HashMap<String, Stat> results = new HashMap<>();

            boolean toAdd = false;

            for (Iterator<String[]> iterator = annotationrows.iterator(); iterator.hasNext();) {
                String[] next = iterator.next();
                if (next[0].contains("SV")) {
                    for (int i = 0; i < next.length; i++) {
                        headerHelper.put(next[i], i);
                    }
                } else if (next[headerHelper.get("AnnotSV type")].equals("full")) {
                    rowbyID.put(next[6], next);
                }
            }

            for (Iterator<String[]> iterator = rows.iterator(); iterator.hasNext();) {
                String[] row = iterator.next();
                if (row[0].startsWith("#CHROM")) {
                    try {
                        for (int i = 9; i < row.length; i++) {
                            sampleList.add(row[i]);
                            results.put(row[i], new Stat());
                        }

                    } catch (Exception e) {
                        System.out.println("Error in format");
                        e.printStackTrace();
                    }

                } else if (!row[0].contains("#")) {
                    if (rowbyID.containsKey(row[2])) {
                        String ID = row[2];

                        switch (filter) {
                            case "ALL":
                                toAdd = true;
                                break;
                            case "GENIC":
                                toAdd = rowbyID.get(ID)[headerHelper.get("Genic")].equals("true");
                                break;
                            case "EXONIC":
                                toAdd = rowbyID.get(ID)[headerHelper.get("Exonic")].equals("true");
                                break;
                            case "NOVEL":
                                toAdd = rowbyID.get(ID)[headerHelper.get("Novel")].equals("true");
                                break;
                            case "NOVEL_GENIC":
                                toAdd = rowbyID.get(ID)[headerHelper.get("Novel")].equals("true") && rowbyID.get(ID)[headerHelper.get("Genic")].equals("true");
                                break;
                            case "NOVEL_EXONIC":
                                toAdd = rowbyID.get(ID)[headerHelper.get("Novel")].equals("true") && rowbyID.get(ID)[headerHelper.get("Exonic")].equals("true");
                                break;
                            default:
                                break;
                        }
                        if (toAdd) {

                            String type = getValuebyKey(row[7], "SVTYPE");

                            Integer size = Integer.valueOf(getValuebyKey(row[7], "END")) - Integer.valueOf(row[1]);

                            Integer genes = 0;

                            if (!rowbyID.get(ID)[headerHelper.get("Gene name")].isEmpty()) {
                                genes = rowbyID.get(ID)[headerHelper.get("Gene name")].split(",").length;
                            }

                            for (int i = 9; i < row.length; i++) {
                                String[] cur_format = row[i].split(":");
                                Stat stat = results.get(sampleList.get(i - 9));

                                switch (type) {
                                    case "DEL":
                                        switch (cur_format[0]) {
                                            case "1/1":
                                                stat.setHom_count_del(stat.getHom_count_del() + 1);
                                                stat.getDel_hom_sizes().add(size);
                                                stat.setGene_number_hom_del(stat.getGene_number_hom_del() + genes);
                                                results.put(sampleList.get(i - 9), stat);
                                                break;
                                            case "0/1":
                                                stat.setHet_count_del(stat.getHet_count_del() + 1);
                                                stat.getDel_het_sizes().add(size);
                                                stat.setGene_number_het_del(stat.getGene_number_het_del() + genes);
                                                results.put(sampleList.get(i - 9), stat);
                                                break;
                                            default:
                                                stat.setRef_count_del(stat.getRef_count_del() + 1);
                                                results.put(sampleList.get(i - 9), stat);
                                                break;
                                        }
                                        break;
                                    case "INV":
                                        switch (cur_format[0]) {
                                            case "1/1":
                                                stat.setHom_count_inv(stat.getHom_count_inv() + 1);
                                                stat.getInv_hom_sizes().add(size);
                                                stat.setGene_number_hom_inv(stat.getGene_number_hom_inv() + genes);
                                                results.put(sampleList.get(i - 9), stat);
                                                break;
                                            case "0/1":
                                                stat.setHet_count_inv(stat.getHet_count_inv() + 1);
                                                stat.getInv_het_sizes().add(size);
                                                stat.setGene_number_het_inv(stat.getGene_number_het_inv() + genes);
                                                results.put(sampleList.get(i - 9), stat);
                                                break;
                                            default:
                                                stat.setRef_count_inv(stat.getRef_count_inv() + 1);
                                                results.put(sampleList.get(i - 9), stat);
                                                break;
                                        }
                                        break;
                                    case "DUP":
                                        switch (cur_format[0]) {
                                            case "1/1":
                                                stat.setHom_count_dup(stat.getHom_count_dup() + 1);
                                                stat.getDup_hom_sizes().add(size);
                                                stat.setGene_number_hom_dup(stat.getGene_number_hom_dup() + genes);
                                                results.put(sampleList.get(i - 9), stat);
                                                break;
                                            case "0/1":
                                                stat.setHet_count_dup(stat.getHet_count_dup() + 1);
                                                stat.getDup_het_sizes().add(size);
                                                stat.setGene_number_het_dup(stat.getGene_number_het_dup() + genes);
                                                results.put(sampleList.get(i - 9), stat);
                                                break;
                                            default:
                                                stat.setRef_count_dup(stat.getRef_count_dup() + 1);
                                                results.put(sampleList.get(i - 9), stat);
                                                break;
                                        }
                                        break;
                                }

                            }

                        }
                    }
                }
                toAdd = false;
            }
            for (int i = 0; i < sampleList.size(); i++) {
                System.out.println(sampleList.get(i) + "\t" + results.get(sampleList.get(i)).toString());
            }

        } catch (IOException ex) {
            Logger.getLogger(Helper.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public static void phenotypeAnnotation(String vcfFile, String affectedsFile) {
        /* INPUT: VCF FILE, AFFECTED FILE (SAMPLE AFFECTED) 
        
        Generate report that contains: ID    COUNT_CASE_HET  COUNT_CASE_HOM  COUNT_CONTROL_HET   COUNT_CONTROL_HOM
        
         */

        try {
            CSVReader reader = new CSVReader(new FileReader(vcfFile), '\t');
            List<String[]> rows = reader.readAll();

            CSVReader affectedReader = new CSVReader(new FileReader(affectedsFile), '\t');
            List<String[]> affectStatusSampleList = affectedReader.readAll();

            HashMap<String, Integer> sampleAffectedStatus = new HashMap<>();

            ArrayList<String> sampleList = new ArrayList<>();

            affectStatusSampleList.forEach((curSampleStatus) -> {
                sampleAffectedStatus.put(curSampleStatus[0], Integer.valueOf(curSampleStatus[1]));
            });

            System.out.println("ID\tCASE_HET\tCASE_HOM\tCONTROL_HET\tCONTROL_HOM");

            for (String[] row : rows) {
                if (row.length > 1) {
                    if (row[0].contains("#CHROM")) {
                        for (int i = 9; i < row.length; i++) {
                            sampleList.add(row[i]);
                        }
                    } else {
                        int count_case_het = 0;
                        int count_case_hom = 0;
                        int control_case_het = 0;
                        int control_case_hom = 0;
                        for (int i = 9; i < row.length; i++) {
                            if (sampleAffectedStatus.getOrDefault(sampleList.get(i - 9), 1).compareTo(2) == 0) {
                                if (row[i].contains("0/1:")) {
                                    count_case_het = count_case_het + 1;
                                } else if (row[i].contains("1/1:")) {
                                    count_case_hom = count_case_hom + 1;
                                }
                            } else {
                                if (row[i].contains("0/1:")) {
                                    control_case_het = control_case_het + 1;
                                } else if (row[i].contains("1/1:")) {
                                    control_case_hom = control_case_hom + 1;
                                }
                            }

                        }
                        System.out.println(row[2] + "\t" + count_case_het + "\t" + count_case_hom + "\t" + control_case_het + "\t" + control_case_hom);

                    }
                }
            }

        } catch (IOException ex) {
            Logger.getLogger(Helper.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    public static void generateReportsbyPED(String filename, String pedFile, String mappingTSV) {
        try {

            CSVReader reader = new CSVReader(new FileReader(filename), '\t');
            List<String[]> rows = reader.readAll();

            CSVReader pedReader = new CSVReader(new FileReader(pedFile), '\t');
            List<String[]> pedRows = pedReader.readAll();

            CSVReader mappingReader = new CSVReader(new FileReader(mappingTSV), '\t');
            List<String[]> mappingRows = mappingReader.readAll();

            HashMap<String, Integer> header = new HashMap<>();
            HashMap<String, String> mapping = new HashMap<>();

            List<String[]> affecteds = new ArrayList<>();

            for (int i = 0; i < pedRows.size(); i++) {
                if (pedRows.get(i)[5].equals("2")) {
                    String Family = pedRows.get(i)[0];
                    String kidSample = pedRows.get(i)[1];
                    String momSample = pedRows.get(i)[2];
                    String dadSample = pedRows.get(i)[3];

                    affecteds.add(new String[]{Family, kidSample, momSample, dadSample});
                }
            }

            for (int i = 0; i < mappingRows.size(); i++) {
                mapping.put(mappingRows.get(i)[1], mappingRows.get(i)[0]);
            }
            System.out.println("ID\tDENOVO\tHOMOZYGOUS");
            for (Iterator<String[]> iterator = rows.iterator(); iterator.hasNext();) {
                String[] row = iterator.next();
                if (row[0].contains("CHROM")) {
                    try {
                        for (int i = 0; i < row.length; i++) {
                            header.put(row[i], i);
                        }
                    } catch (Exception e) {
                        System.out.println("Error in format");
                        e.printStackTrace();
                    }

                } else if (!row[0].contains("#")) {
                    String denovo = "";
                    String homozygous = "";
                    for (int i = 0; i < affecteds.size(); i++) {

                        Integer kidIndex = header.getOrDefault(mapping.get(affecteds.get(i)[1]), 0);
                        Integer momIndex = header.getOrDefault(mapping.get(affecteds.get(i)[2]), 0);
                        Integer dadIndex = header.getOrDefault(mapping.get(affecteds.get(i)[3]), 0);
                        if (kidIndex != 0 && momIndex != 0 && dadIndex != 0) {
                            if (row[kidIndex].contains("0/1:") && row[momIndex].contains("0/0:") && row[dadIndex].contains("0/0:")) {
                                denovo = denovo.concat(affecteds.get(i)[0] + ":" + affecteds.get(i)[1] + "|");
                            } else if (row[kidIndex].contains("1/1:") && row[momIndex].contains("0/1:") && row[dadIndex].contains("0/1:")) {
                                homozygous = homozygous.concat(affecteds.get(i)[0] + ":" + affecteds.get(i)[1] + "|");
                            }
                        }

                    }
                    if (!denovo.isEmpty() || !homozygous.isEmpty()) {
                        System.out.println(row[2] + "\t" + denovo + "\t" + homozygous);
                    } else {
                        System.out.println(row[2] + "\t" + "\t");
                    }

                }

            }

        } catch (IOException ex) {
            Logger.getLogger(Helper.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public static String getValuebyKey(String infoRow, String key) {
        String[] infoValues = infoRow.split(";");
        for (int i = 0; i < infoValues.length; i++) {
            if (infoValues[i].contains(key)) {
                if (infoValues[i].split("=")[0].equals(key)) {
                    return infoValues[i].split("=")[1];
                }
            }
        }
        return "0";
    }

    static void mergePhenotypes(String mapping, String pheno) {
        try {

            List<String[]> ultimateTable = new ArrayList<>();

            ArrayList<String> headerLine = new ArrayList<>();

            List<HashMap> phenoList = new ArrayList<>();

            CSVReader mappingReader = new CSVReader(new FileReader(mapping), '\t');
            List<String[]> mappings = mappingReader.readAll();

            Scanner s = new Scanner(new File("filepath"));
            ArrayList<String> phenolist = new ArrayList<>();
            while (s.hasNext()) {
                phenolist.add(s.next());
            }
            s.close();

            for (int i = 0; i < phenolist.size(); i++) {
                CSVReader phenoReader = new CSVReader(new FileReader(phenolist.get(i)), '\t');
                List<String[]> phenotype = phenoReader.readAll();
                HashMap<String, String[]> phenoMap = new HashMap<>();
                for (int j = 0; j < phenotype.size(); j++) {
                    if (phenotype.get(j)[0].equals("DummyID")) {
                        headerLine.addAll(Arrays.asList(phenotype.get(j)));
                    } else {
                        String[] currow = phenoMap.get(0);
                        phenoMap.put(currow[0], new String[]{currow[2], currow[3]});
                    }

                }
                phenoList.add(phenoMap);
            }

            ArrayList<String> Ids = new ArrayList<>();
            for (int a = 0; a < mappings.size(); a++) {
                String ID = mappings.get(a)[5];
                ArrayList<String> ultimateRow = new ArrayList<>();
                for (int i = 0; i < phenoList.size(); i++) {
                    String[] curValues = (String[]) phenoList.get(i).get(ID);
                    ultimateRow.addAll(Arrays.asList(mappings.get(a)));
                    ultimateRow.addAll(Arrays.asList(curValues));
                }
                ultimateTable.add((String[]) ultimateRow.toArray());
                ultimateRow.clear();
            }

        } catch (FileNotFoundException ex) {
            Logger.getLogger(Helper.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(Helper.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    static void annotateFilewithSubPopFrequencies(String vcfFile, String populationStructureClusters) {
        try {

            Reader reader = Files.newBufferedReader(Paths.get(vcfFile));
            CSVReader csvReader = new CSVReader(reader, '\t');

            //Generating pop structure mapping 
            CSVReader popLabels = new CSVReader(new FileReader(populationStructureClusters), '\t');
            List<String[]> popRows = popLabels.readAll();
            HashMap<String, String> popMapping = new HashMap<>();

            List<String[]> output = new ArrayList<>();

            for (String[] next : popRows) {
                popMapping.put(next[0], next[1]);
            }

            ArrayList<String> sampleList = new ArrayList<>();

            String[] row;

            String replace = vcfFile.replace(".vcf", ".subpop.vcf");
            Writer writer = Files.newBufferedWriter(Paths.get(replace));
            CSVWriter csvwriter = new CSVWriter(writer, '\t', CSVWriter.NO_QUOTE_CHARACTER);

            while ((row = csvReader.readNext()) != null) {
                if (row[0].startsWith("#CHROM")) {
                    csvwriter.writeNext(row);
                    try {
                        for (int i = 9; i < row.length; i++) {
                            sampleList.add(row[i]);
                        }
                    } catch (Exception e) {
                        System.out.println("Error in format" + e);
                    }

                } else if (!row[0].contains("#")) {
                    //Introducting QGP Labels
                    Integer qgp_adm_ac = 0;
                    Integer qgp_adm_an = 0;
                    Integer qgp_adm_homo = 0;
                    double qgp_adm_af;
                    Integer qgp_afr_ac = 0;
                    Integer qgp_afr_an = 0;
                    Integer qgp_afr_homo = 0;
                    double qgp_afr_af;
                    Integer qgp_gar_ac = 0;
                    Integer qgp_gar_an = 0;
                    Integer qgp_gar_homo = 0;
                    double qgp_gar_af;
                    Integer qgp_par_ac = 0;
                    Integer qgp_par_an = 0;
                    Integer qgp_par_homo = 0;
                    double qgp_par_af;
                    Integer qgp_sas_ac = 0;
                    Integer qgp_sas_an = 0;
                    Integer qgp_sas_homo = 0;
                    double qgp_sas_af;
                    Integer qgp_wep_ac = 0;
                    Integer qgp_wep_an = 0;
                    Integer qgp_wep_homo = 0;
                    double qgp_wep_af;

                    for (int i = 9; i < row.length; i++) {
                        String[] cur_format = row[i].split(":");

                        String popColor = popMapping.getOrDefault(sampleList.get(i - 9), "notfounbd");
                        switch (popColor) {
                            case "Coral":
                                switch (cur_format[0]) {
                                    case "1/1":
                                        qgp_par_ac = qgp_par_ac + 2;
                                        qgp_par_an = qgp_par_an + 2;
                                        qgp_par_homo = qgp_par_homo + 1;
                                        break;
                                    case "0/1":
                                        qgp_par_ac = qgp_par_ac + 1;
                                        qgp_par_an = qgp_par_an + 2;
                                        break;
                                    case "0/0":
                                        qgp_par_an = qgp_par_an + 2;
                                        break;
                                    default:
                                        break;
                                }
                                break;
                            case "Blueberry":
                                switch (cur_format[0]) {
                                    case "1/1":
                                        qgp_gar_ac = qgp_gar_ac + 2;
                                        qgp_gar_an = qgp_gar_an + 2;
                                        qgp_gar_homo = qgp_gar_homo + 1;
                                        break;
                                    case "0/1":
                                        qgp_gar_ac = qgp_gar_ac + 1;
                                        qgp_gar_an = qgp_gar_an + 2;
                                        break;
                                    case "0/0":
                                        qgp_gar_an = qgp_gar_an + 2;
                                        break;
                                    default:
                                        break;
                                }
                                break;
                            case "Strawberry":
                                switch (cur_format[0]) {
                                    case "1/1":
                                        qgp_wep_ac = qgp_wep_ac + 2;
                                        qgp_wep_an = qgp_wep_an + 2;
                                        qgp_wep_homo = qgp_wep_homo + 1;
                                        break;
                                    case "0/1":
                                        qgp_wep_ac = qgp_wep_ac + 1;
                                        qgp_wep_an = qgp_wep_an + 2;
                                        break;
                                    case "0/0":
                                        qgp_wep_an = qgp_wep_an + 2;
                                        break;
                                    default:
                                        break;
                                }
                                break;
                            case "Lemon":
                                switch (cur_format[0]) {
                                    case "1/1":
                                        qgp_sas_ac = qgp_sas_ac + 2;
                                        qgp_sas_an = qgp_sas_an + 2;
                                        qgp_sas_homo = qgp_sas_homo + 1;
                                        break;
                                    case "0/1":
                                        qgp_sas_ac = qgp_sas_ac + 1;
                                        qgp_sas_an = qgp_sas_an + 2;
                                        break;
                                    case "0/0":
                                        qgp_sas_an = qgp_sas_an + 2;
                                        break;
                                    default:
                                        break;
                                }
                                break;
                            case "Orange":
                                switch (cur_format[0]) {
                                    case "1/1":
                                        qgp_afr_ac = qgp_afr_ac + 2;
                                        qgp_afr_an = qgp_afr_an + 2;
                                        qgp_afr_homo = qgp_afr_homo + 1;
                                        break;
                                    case "0/1":
                                        qgp_afr_ac = qgp_afr_ac + 1;
                                        qgp_afr_an = qgp_afr_an + 2;
                                        break;
                                    case "0/0":
                                        qgp_afr_an = qgp_afr_an + 2;
                                        break;
                                    default:
                                        break;
                                }
                                break;
                            case "Mushroom":
                                switch (cur_format[0]) {
                                    case "1/1":
                                        qgp_adm_ac = qgp_adm_ac + 2;
                                        qgp_adm_an = qgp_adm_an + 2;
                                        qgp_adm_homo = qgp_adm_homo + 1;
                                        break;
                                    case "0/1":
                                        qgp_adm_ac = qgp_adm_ac + 1;
                                        qgp_adm_an = qgp_adm_an + 2;
                                        break;
                                    case "0/0":
                                        qgp_adm_an = qgp_adm_an + 2;
                                        break;
                                    default:
                                        break;
                                }
                                break;

                        }

                    }
                    qgp_par_af = ((double) qgp_par_ac / ((double) (qgp_par_an)));
                    qgp_gar_af = ((double) qgp_gar_ac / ((double) (qgp_gar_an)));
                    qgp_wep_af = ((double) qgp_wep_ac / ((double) (qgp_wep_an)));
                    qgp_sas_af = ((double) qgp_sas_ac / ((double) (qgp_sas_an)));
                    qgp_afr_af = ((double) qgp_afr_ac / ((double) (qgp_afr_an)));
                    qgp_adm_af = ((double) qgp_adm_ac / ((double) (qgp_adm_an)));

                    row[7] = row[7] + ";QGP_ADM_AC=" + qgp_adm_ac + ";"
                            + "QGP_ADM_AF=" + qgp_adm_af + ";"
                            + "QGP_ADM_AN=" + qgp_adm_an + ";"
                            + "QGP_ADM_HOMO=" + qgp_adm_homo + ";"
                            + "QGP_AFR_AC=" + qgp_afr_ac + ";"
                            + "QGP_AFR_AF=" + qgp_afr_af + ";"
                            + "QGP_AFR_AN=" + qgp_afr_an + ";"
                            + "QGP_AFR_HOMO=" + qgp_afr_homo + ";"
                            + "QGP_GAR_AC=" + qgp_gar_ac + ";"
                            + "QGP_GAR_AF=" + qgp_gar_af + ";"
                            + "QGP_GAR_AN=" + qgp_gar_an + ";"
                            + "QGP_GAR_HOMO=" + qgp_gar_homo + ";"
                            + "QGP_PAR_AC=" + qgp_par_ac + ";"
                            + "QGP_PAR_AF=" + qgp_par_af + ";"
                            + "QGP_PAR_AN=" + qgp_par_an + ";"
                            + "QGP_PAR_HOMO=" + qgp_par_homo + ";"
                            + "QGP_SAS_AC=" + qgp_sas_ac + ";"
                            + "QGP_SAS_AF=" + qgp_sas_af + ";"
                            + "QGP_SAS_AN=" + qgp_sas_an + ";"
                            + "QGP_SAS_HOMO=" + qgp_sas_homo + ";"
                            + "QGP_WEP_AC=" + qgp_wep_ac + ";"
                            + "QGP_WEP_AF=" + qgp_wep_af + ";"
                            + "QGP_WEP_AN=" + qgp_wep_an + ";"
                            + "QGP_WEP_HOMO=" + qgp_wep_homo + ";";
                    csvwriter.writeNext(row);

                }
            }
            Thread.sleep(2000);
            csvwriter.close();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(Helper.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(Helper.class.getName()).log(Level.SEVERE, null, ex);
        } catch (InterruptedException ex) {
            Logger.getLogger(Helper.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    static void ranker(String tsvFile, String outputFile) throws FileNotFoundException, IOException, InterruptedException {
        CSVReader popLabels = new CSVReader(new FileReader(tsvFile), ',');
        List<String[]> tsvPhenoFile = popLabels.readAll();
        String[] header = tsvPhenoFile.get(0);
        tsvPhenoFile.remove(0);

        Collections.sort(tsvPhenoFile, new Comparator<String[]>() {
            @Override
            public int compare(String[] o1, String[] o2) {
                if (o1.length == 0) {
                    return o2.length == 0 ? 0 : -1;
                }
                if (o2.length == 0) {
                    return 1;
                }
                return Double.valueOf(o2[2]).compareTo(Double.valueOf(o1[2]));
            }
        });

        int len = tsvPhenoFile.size();
        Map<String, int[]> map = new TreeMap<>((a, b) -> Double.valueOf(b).compareTo(Double.valueOf(a)));
        for (int i = 0; i < len; i++) {
            String key = tsvPhenoFile.get(i)[2];
            int[] arr;
            if ((arr = map.get(key)) == null) {
                arr = new int[]{i, 0};
            }
            arr[1]++;
            map.put(key, arr);
        }

        int[][] groups = map.values().toArray(new int[map.size()][]);
        standardRanking(len, groups, tsvPhenoFile, tsvFile, header, outputFile);
    }

    private static void standardRanking(int len, int[][] groups, List<String[]> list, String tsvFile, String[] header, String outputFile) throws IOException, InterruptedException {
        System.out.println("\nStandard ranking");

        String[] rankArrayHeader = {String.valueOf(tsvFile.replace(" ", "") + "_Rank")};
        String[] headerRow = (String[]) ArrayUtils.addAll(header, rankArrayHeader);

        String replace = tsvFile.replace(".tsv", ".ranked.tsv");
        Writer writer = Files.newBufferedWriter(Paths.get(outputFile));
        CSVWriter csvwriter = new CSVWriter(writer, '\t', CSVWriter.NO_QUOTE_CHARACTER);

        csvwriter.writeNext(headerRow);

        ArrayList<String[]> output = new ArrayList<>();
        for (int i = 0, rank = 0, group = 0; i < len; i++) {
            if (group < groups.length && i == groups[group][0]) {
                rank = i + 1;
                group++;
            }
            String[] rankArray = {String.valueOf(rank)};
            String[] newRow = (String[]) ArrayUtils.addAll(list.get(i), rankArray);
            csvwriter.writeNext(newRow);
        }
        Thread.sleep(2000);
        csvwriter.close();
    }

    static void dummyVCF(String regions, String outputFile) throws FileNotFoundException, IOException, InterruptedException {

        CSVReader reader = new CSVReader(new FileReader(regions), '\t');
        List<String[]> dummyRegions = reader.readAll();

        List<String[]> dummyVCF = new ArrayList<>();

        dummyVCF.add(new String[]{"##fileformat=VCFv4.1"});
        dummyVCF.add(new String[]{"##source=SURVIVOR"});
        dummyVCF.add(new String[]{"##fileDate=20190324"});
        dummyVCF.add(new String[]{"##contig=<ID=1,length=249250621>"});
        dummyVCF.add(new String[]{"##contig=<ID=10,length=135534747>"});
        dummyVCF.add(new String[]{"##contig=<ID=11,length=135006516>"});
        dummyVCF.add(new String[]{"##contig=<ID=12,length=133851895>"});
        dummyVCF.add(new String[]{"##contig=<ID=13,length=115169878>"});
        dummyVCF.add(new String[]{"##contig=<ID=14,length=107349540>"});
        dummyVCF.add(new String[]{"##contig=<ID=15,length=102531392>"});
        dummyVCF.add(new String[]{"##contig=<ID=16,length=90354753>"});
        dummyVCF.add(new String[]{"##contig=<ID=17,length=81195210>"});
        dummyVCF.add(new String[]{"##contig=<ID=18,length=78077248>"});
        dummyVCF.add(new String[]{"##contig=<ID=19,length=59128983>"});
        dummyVCF.add(new String[]{"##contig=<ID=2,length=243199373>"});
        dummyVCF.add(new String[]{"##contig=<ID=20,length=63025520>"});
        dummyVCF.add(new String[]{"##contig=<ID=21,length=48129895>"});
        dummyVCF.add(new String[]{"##contig=<ID=22,length=51304566>"});
        dummyVCF.add(new String[]{"##contig=<ID=3,length=198022430>"});
        dummyVCF.add(new String[]{"##contig=<ID=4,length=191154276>"});
        dummyVCF.add(new String[]{"##contig=<ID=5,length=180915260>"});
        dummyVCF.add(new String[]{"##contig=<ID=6,length=171115067>"});
        dummyVCF.add(new String[]{"##contig=<ID=7,length=159138663>"});
        dummyVCF.add(new String[]{"##contig=<ID=8,length=146364022>"});
        dummyVCF.add(new String[]{"##contig=<ID=9,length=141213431>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000191.1,length=106433>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000192.1,length=547496>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000193.1,length=189789>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000194.1,length=191469>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000195.1,length=182896>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000196.1,length=38914>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000197.1,length=37175>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000198.1,length=90085>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000199.1,length=169874>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000200.1,length=187035>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000201.1,length=36148>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000202.1,length=40103>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000203.1,length=37498>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000204.1,length=81310>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000205.1,length=174588>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000206.1,length=41001>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000207.1,length=4262>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000208.1,length=92689>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000209.1,length=159169>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000210.1,length=27682>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000211.1,length=166566>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000212.1,length=186858>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000213.1,length=164239>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000214.1,length=137718>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000215.1,length=172545>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000216.1,length=172294>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000217.1,length=172149>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000218.1,length=161147>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000219.1,length=179198>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000220.1,length=161802>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000221.1,length=155397>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000222.1,length=186861>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000223.1,length=180455>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000224.1,length=179693>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000225.1,length=211173>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000226.1,length=15008>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000227.1,length=128374>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000228.1,length=129120>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000229.1,length=19913>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000230.1,length=43691>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000231.1,length=27386>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000232.1,length=40652>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000233.1,length=45941>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000234.1,length=40531>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000235.1,length=34474>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000236.1,length=41934>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000237.1,length=45867>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000238.1,length=39939>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000239.1,length=33824>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000240.1,length=41933>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000241.1,length=42152>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000242.1,length=43523>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000243.1,length=43341>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000244.1,length=39929>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000245.1,length=36651>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000246.1,length=38154>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000247.1,length=36422>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000248.1,length=39786>"});
        dummyVCF.add(new String[]{"##contig=<ID=GL000249.1,length=38502>"});
        dummyVCF.add(new String[]{"##contig=<ID=MT,length=16569>"});
        dummyVCF.add(new String[]{"##contig=<ID=NC_007605,length=171823>"});
        dummyVCF.add(new String[]{"##contig=<ID=X,length=155270560>"});
        dummyVCF.add(new String[]{"##contig=<ID=Y,length=59373566>"});
        dummyVCF.add(new String[]{"##contig=<ID=hs37d5,length=35477943>"});
        dummyVCF.add(new String[]{"##ALT=<ID=DEL,Description=\"Deletion\">"});
        dummyVCF.add(new String[]{"##ALT=<ID=DUP,Description=\"Duplication\">"});
        dummyVCF.add(new String[]{"##ALT=<ID=INV,Description=\"Inversion\">"});
        dummyVCF.add(new String[]{"##ALT=<ID=BND,Description=\"Translocation\">"});
        dummyVCF.add(new String[]{"##ALT=<ID=INS,Description=\"Insertion\">"});
        dummyVCF.add(new String[]{"##INFO=<ID=CIEND,Number=1,Type=String,Description=\"PE confidence interval around END\">"});
        dummyVCF.add(new String[]{"##INFO=<ID=CIPOS,Number=1,Type=String,Description=\"PE confidence interval around POS\">"});
        dummyVCF.add(new String[]{"##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">"});
        dummyVCF.add(new String[]{"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">"});
        dummyVCF.add(new String[]{"##INFO=<ID=MAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of paired-ends\">"});
        dummyVCF.add(new String[]{"##INFO=<ID=RE,Number=1,Type=Integer,Description=\"read support\">"});
        dummyVCF.add(new String[]{"##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">"});
        dummyVCF.add(new String[]{"##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">"});
        dummyVCF.add(new String[]{"##INFO=<ID=AVGLEN,Number=1,Type=Float,Description=\"Length of the SV\">"});
        dummyVCF.add(new String[]{"##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Method for generating this merged VCF file.\">"});
        dummyVCF.add(new String[]{"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV.\">"});
        dummyVCF.add(new String[]{"##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description=\"Vector of supporting samples.\">"});
        dummyVCF.add(new String[]{"##INFO=<ID=SUPP,Number=1,Type=String,Description=\"Number of samples supporting the variant\">"});
        dummyVCF.add(new String[]{"##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Indicating the direction of the reads with respect to the type and breakpoint.\">"});
        dummyVCF.add(new String[]{"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"});
        dummyVCF.add(new String[]{"##FORMAT=<ID=LN,Number=1,Type=Integer,Description=\"predicted length\">"});
        dummyVCF.add(new String[]{"##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# supporting reference,variant reads in that order\">"});
        dummyVCF.add(new String[]{"##FORMAT=<ID=ST,Number=1,Type=String,Description=\"Strand of SVs\">"});
        dummyVCF.add(new String[]{"##FORMAT=<ID=TY,Number=1,Type=String,Description=\"Types\">"});
        dummyVCF.add(new String[]{"##FORMAT=<ID=CO,Number=1,Type=String,Description=\"Coordinates\">"});
        dummyVCF.add(new String[]{"##FORMAT=<ID=PSV,Number=1,Type=String,Description=\"Previous support vector\">"});
        //0 - CHROM 1 - POS 2 - END 3 - SIZE 4 - TYPE
        dummyVCF.add(new String[]{"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"});

        for (int i = 0; i < dummyRegions.size(); i++) {

            String INFO = "SUPP=3;SUPP_VEC=111;AVGLEN=" + dummyRegions.get(i)[3]
                    + ";SVTYPE=" + dummyRegions.get(i)[4] + ";SVMETHOD=SURVIVORv2;CHR2=1;END=" + dummyRegions.get(i)[2] + ";CIPOS=-1000,1000;CIEND=-1000,1000;STRANDS=+-";

            String[] cnvregion = {dummyRegions.get(i)[0], dummyRegions.get(i)[1], String.valueOf(i), "N", "<" + dummyRegions.get(i)[4] + ">", ".", "PASS", INFO, "GT:PSV:LN:DR:ST:TY:CO"};
            dummyVCF.add(cnvregion);
        }

        System.out.println(dummyVCF.size());

        CSVWriter writer = new CSVWriter(new FileWriter(outputFile), '\t', CSVWriter.NO_QUOTE_CHARACTER, CSVWriter.NO_QUOTE_CHARACTER, CSVWriter.DEFAULT_LINE_END);
        writer.writeAll(dummyVCF);
        Thread.sleep(2000);
        writer.close();

    }

    static void phenoCalculator(String vcfFile, int vcfIdColumn, String phenoFile, int phenoIdcolumn, String outputFile) throws IOException, InterruptedException {
        CSVReader cnvReader = new CSVReader(new FileReader(vcfFile), '\t');
        List<String[]> cnvs = cnvReader.readAll();
        cnvReader.close();

        CSVReader phenoReader = new CSVReader(new FileReader(phenoFile), '\t');
        List<String[]> phenos = phenoReader.readAll();
        phenoReader.close();

        List<Integer> phenoIds = new ArrayList<>();
        List<String> sampleList = new ArrayList<>();

        HashMap<String, String[]> phenoTypeMap = new HashMap<>();

        for (int i = 0; i < phenos.size(); i++) {
            phenoTypeMap.put(phenos.get(i)[phenoIdcolumn], phenos.get(i));
        }

        for (int i = 0; i < cnvs.size(); i++) {
            String[] currentRow = cnvs.get(i);

            if (currentRow[0].startsWith("#CHROM")) {
                for (int j = 9; j < currentRow.length; j++) {
                    sampleList.add(currentRow[j]);
                }
            }
        }

        for (int i = 0; i < cnvs.size(); i++) {
            ArrayList<Double> refMeanList = new ArrayList<>();
            ArrayList<Double> hetMeanList = new ArrayList<>();
            ArrayList<Double> homMeanList = new ArrayList<>();

            String[] currentRow = cnvs.get(i);
            if (!currentRow[0].startsWith("#")) {
                for (int j = 9; j < currentRow.length; j++) {
                    if (currentRow[j].contains("0/0")) {
                        refMeanList.add(Double.valueOf(currentRow[j]));
                    } else if (currentRow[j].contains("0/1")) {
                        hetMeanList.add(Double.valueOf(currentRow[j]));
                    } else if (currentRow[j].contains("1/1")) {
                        homMeanList.add(Double.valueOf(currentRow[j]));
                    }
                }

            }
        }

        List<String[]> cnvMeanReport = new ArrayList<>();

        CSVWriter writer = new CSVWriter(new FileWriter(outputFile), '\t', CSVWriter.NO_QUOTE_CHARACTER);
        writer.writeAll(cnvMeanReport);
        Thread.sleep(2000);
        writer.close();

    }

    public static double mean(double[] m) {
        double sum = 0;
        for (int i = 0; i < m.length; i++) {
            sum += m[i];
        }
        return sum / m.length;
    }

    public static void svInspector(String bamFile, String outFolder, String outName, String regions, String referencePath, String bcftoolsPath, String svtyperPath, String dupholdPath) {
        try {
            String name = bamFile.substring(bamFile.lastIndexOf("/") + 1).replace(".bwa.bam", "");
            String fullPath = outFolder + name;
            Helper.dummyVCF(regions, fullPath + ".dummy.vcf");

            ProcessBuilder processBuilder = new ProcessBuilder();

            processBuilder.command(svtyperPath, "-i", fullPath + ".dummy.vcf", "-B", bamFile);
            processBuilder.redirectOutput(new File(fullPath + ".gt.vcf"));
            Process process = processBuilder.start();
            addOutputListener(process);
            processBuilder.command(dupholdPath, "-v", fullPath + ".gt.vcf", "-b", bamFile, "-f", referencePath, "-t", "4", "-o", fullPath + "gt.duphold.vcf");
            processBuilder.redirectOutput();
            Process dupholdprocess = processBuilder.start();
            addOutputListener(dupholdprocess);


            /*
            //./duphold -v $2/$name.gt.vcf -b $line -f $ref_data -t 4 -o $2/$name.gt.duphold.vcf
            processBuilder.command(dupholdPath, "-v", fullPath + "gt.vcf", "-b", bamFile, "-f", referencePath, "-t", "4", "-o", fullPath + "gt.duphold.vcf");
            //bcftools query -f '[%SAMPLE]\t%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/END\t[%GT]\t[%DHFC]\t[%DHBFC]\t[%DHFFC]\n' $2/$name.gt.duphold.vcf >> $2/$3
            processBuilder.command(bcftoolsPath,
                    "query -f",
                    "'[%SAMPLE]\\t%CHROM\\t%POS\\t%INFO/END\\t%INFO/SVTYPE\\t%INFO/END\\t[%GT]\\t[%DHFC]\\t[%DHBFC]\\t[%DHFFC]\\n'",
                    fullPath + "gt.duphold.vcf",
                    ">>", outFolder + outName);
            //rm -rf $2/$name.dummy.vcf;
            processBuilder.command("rm", "-rf", "", fullPath + "dummy.vcf");
            //rm -rf $2/$name.gt.vcf;
            processBuilder.command("rm", "-rf", "", fullPath + "gt.vcf");*/
        } catch (IOException | InterruptedException ex) {
            Logger.getLogger(Helper.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    private static void addOutputListener(Process process) throws InterruptedException, IOException {
        StringBuilder output = new StringBuilder();

        BufferedReader reader = new BufferedReader(
                new InputStreamReader(process.getErrorStream()));

        String line;
        while ((line = reader.readLine()) != null) {
            output.append(line).append("\n");
        }

        int exitVal = process.waitFor();
        if (exitVal == 0) {
            System.out.println("Success!");
            System.out.println(output);
        } else {
            System.out.println("Failed!");
            System.out.println(output);
        }
    }

    static void clusterVCF(String vcfFile) {
        Reader reader = null;
        try {
            CSVReader cnvReader = new CSVReader(new FileReader(vcfFile), '\t');
            List<String[]> cnvs = cnvReader.readAll();
            cnvReader.close();
            ArrayList<String> sampleList = new ArrayList<>();

            String[] row;

            String replace = vcfFile.replace(".vcf", ".mcnv.vcf");
            Writer writer = Files.newBufferedWriter(Paths.get(replace));
            CSVWriter csvwriter = new CSVWriter(writer, '\t', CSVWriter.NO_QUOTE_CHARACTER);

            for (int i = 0; i < cnvs.size(); i++) {
                row = cnvs.get(i);
                if (row[0].startsWith("#CHROM")) {
                    csvwriter.writeNext(row);
                    try {
                        for (int j = 9; j < row.length; j++) {
                            sampleList.add(row[i]);
                        }
                    } catch (Exception e) {
                        System.out.println("Error in format");
                        e.printStackTrace();
                    }

                } else if (!row[0].contains("#")) {

                }
            }

        } catch (IOException ex) {
            Logger.getLogger(Helper.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                reader.close();
            } catch (IOException ex) {
                Logger.getLogger(Helper.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

    }
    

    static void annotateVCFwithPhenos(String vcfFile, String phenofile) {

        try {
            Reader reader = Files.newBufferedReader(Paths.get(vcfFile));
            CSVReader csvReader = new CSVReader(reader, '\t');

            CSVReader phenoReader = new CSVReader(new FileReader(phenofile), '\t');
            List<String[]> phenoRows = phenoReader.readAll();
            ArrayList<Phenotype> phenotypes = new ArrayList<>();
            for (int i = 0; i < phenoRows.size(); i++) {
                if (phenoRows.get(i)[0].contains("Dummy")) {
                    for (int j = 0; j < phenoRows.get(i).length / 6; j++) {
                        phenotypes.add(new Phenotype(phenoRows.get(i)[1 + (j * 6)], new HashMap<>(), new HashMap<>(), new HashMap<>()));
                    }
                } else {
                    for (int j = 0; j < phenotypes.size(); j++) {
                        phenotypes.get(j).getPheno_raw().put(phenoRows.get(i)[0], Double.valueOf(phenoRows.get(i)[1 + (j * 6)]));
                        phenotypes.get(j).getPheno_filtered().put(phenoRows.get(i)[0], Double.valueOf(phenoRows.get(i)[2 + (j * 6)]));
                        phenotypes.get(j).getPheno_Rank().put(phenoRows.get(i)[0], Double.valueOf(phenoRows.get(i)[5 + (j * 6)]));

                    }
                }

            }

            ArrayList<String> sampleList = new ArrayList<>();

            String[] row;
            ArrayList<String> outputRow = new ArrayList<>();

            String replace = vcfFile.replace(".vcf", ".pheno.tsv");
            Writer writer = Files.newBufferedWriter(Paths.get(replace));
            CSVWriter csvwriter = new CSVWriter(writer, '\t', CSVWriter.NO_QUOTE_CHARACTER);
            ArrayList refSamples = new ArrayList();
            ArrayList hetSamples = new ArrayList();
            ArrayList homSamples = new ArrayList();

            while ((row = csvReader.readNext()) != null) {
                refSamples.clear();
                hetSamples.clear();
                homSamples.clear();
                outputRow.clear();
                if (row[0].startsWith("#CHROM")) {
                    try {
                        for (int i = 9; i < row.length; i++) {
                            sampleList.add(row[i]);
                        }
                    } catch (Exception e) {
                        System.out.println("Error in format" + e);
                    }

                    outputRow.add("SV_ID");

                    for (Phenotype phenotype : phenotypes) {
                        outputRow.add(phenotype.getName() + "REF_Count");
                        outputRow.add(phenotype.getName() + "HET_Count");
                        outputRow.add(phenotype.getName() + "HOM_Count");
                        outputRow.add(phenotype.getName() + "REF_Raw_MN");
                        outputRow.add(phenotype.getName() + "HET_Raw_MN");
                        outputRow.add(phenotype.getName() + "HOM_Raw_MN");
                        outputRow.add(phenotype.getName() + "REF_Filtered_MN");
                        outputRow.add(phenotype.getName() + "HET_Filtered_MN");
                        outputRow.add(phenotype.getName() + "HOM_Filtered_MN");
                        outputRow.add(phenotype.getName() + "REF_Rank_MN");
                        outputRow.add(phenotype.getName() + "HET_Rank_MN");
                        outputRow.add(phenotype.getName() + "HOM_Rank_MN");
                    }

                } else if (!row[0].contains("#")) {
                    for (int i = 9; i < row.length; i++) {
                        if (row[i].contains("0/1:")) {
                            hetSamples.add(sampleList.get(i - 9));
                        } else if (row[i].contains("1/1:")) {
                            homSamples.add(sampleList.get(i - 9));
                        } else {
                            refSamples.add(sampleList.get(i - 9));
                        }
                    }

                    outputRow.add(row[2]);

                    for (int i = 0; i < phenotypes.size(); i++) {
                        outputRow.add(String.valueOf(refSamples.size()));
                        outputRow.add(String.valueOf(hetSamples.size()));
                        outputRow.add(String.valueOf(homSamples.size()));
                        outputRow.add(String.valueOf(calculateAverage(refSamples, phenotypes.get(i).getPheno_raw())));
                        outputRow.add(String.valueOf(calculateAverage(hetSamples, phenotypes.get(i).getPheno_raw())));
                        outputRow.add(String.valueOf(calculateAverage(homSamples, phenotypes.get(i).getPheno_raw())));
                        outputRow.add(String.valueOf(calculateAverage(refSamples, phenotypes.get(i).getPheno_filtered())));
                        outputRow.add(String.valueOf(calculateAverage(hetSamples, phenotypes.get(i).getPheno_filtered())));
                        outputRow.add(String.valueOf(calculateAverage(homSamples, phenotypes.get(i).getPheno_filtered())));
                        outputRow.add(String.valueOf(calculateAverage(refSamples, phenotypes.get(i).getPheno_Rank())));
                        outputRow.add(String.valueOf(calculateAverage(hetSamples, phenotypes.get(i).getPheno_Rank())));
                        outputRow.add(String.valueOf(calculateAverage(homSamples, phenotypes.get(i).getPheno_Rank())));
                    }

                }
                String[] output = outputRow.stream().toArray(String[]::new);
                csvwriter.writeNext(output);
            }
            Thread.sleep(2000);
            writer.close();

        } catch (FileNotFoundException ex) {
            Logger.getLogger(Helper.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException | InterruptedException ex) {
            Logger.getLogger(Helper.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    static void adjustGenotype(String vcfFile) {

        Reader reader = null;
        try {
            reader = Files.newBufferedReader(Paths.get(vcfFile));
            CSVReader csvReader = new CSVReader(reader, '\t');

            String replace = vcfFile.replace(".vcf", ".filtered.vcf");
            Writer writer = Files.newBufferedWriter(Paths.get(replace));
            CSVWriter csvwriter = new CSVWriter(writer, '\t', CSVWriter.NO_QUOTE_CHARACTER);

            String[] row;

            while ((row = csvReader.readNext()) != null) {
                if (row[0].startsWith("#CHROM")) {
                    csvwriter.writeNext(row);
                } else {
                    for (int i = 9; i < row.length; i++) {
                        String genotypeQuality = row[i].split(":")[2];
                        if (genotypeQuality.equals("LQ")) {
                            row[i] = row[i].replace("1/1:", "0/0:");
                            row[i] = row[i].replace("0/1:", "0/0:");
                        }
                    }

                    csvwriter.writeNext(row);
                }
            }

        } catch (IOException ex) {
            Logger.getLogger(Helper.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                Thread.sleep(2000);
                reader.close();
            } catch (IOException | InterruptedException ex) {
                Logger.getLogger(Helper.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    static void adjustALT(String vcfFile) {

        Reader reader = null;
        try {
            reader = Files.newBufferedReader(Paths.get(vcfFile));
            CSVReader csvReader = new CSVReader(reader, '\t');

            String replace = vcfFile.replace(".vcf", ".altadjusted.vcf");
            Writer writer = Files.newBufferedWriter(Paths.get(replace));
            CSVWriter csvwriter = new CSVWriter(writer, '\t', CSVWriter.NO_QUOTE_CHARACTER);

            String[] row;

            while ((row = csvReader.readNext()) != null) {
                if (row[0].startsWith("#")) {
                    csvwriter.writeNext(row);
                } else {
                    System.out.println(getValuebyKey(row[7], "SVTYPE"));
                    if (null != getValuebyKey(row[7], "SVTYPE")) {
                        switch (getValuebyKey(row[7], "SVTYPE")) {
                            case "DEL":
                                row[4] = "<DEL>";
                                break;
                            case "DUP":
                                row[4] = "<DUP>";
                                break;
                            case "INV":
                                row[4] = "<INV>";
                                break;
                            case "INS":
                                row[4] = "<INS>";
                                break;
                            default:
                                row[4] = "<TRA>";
                                break;
                        }
                    }

                    csvwriter.writeNext(row);
                }
            }

            Thread.sleep(10000);
            writer.flush();
            writer.close();

        } catch (IOException | InterruptedException ex) {
            Logger.getLogger(Helper.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                reader.close();
                
            } catch (IOException ex) {
                Logger.getLogger(Helper.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    private static double calculateAverage(ArrayList gtSamples, HashMap<String, Double> phenoMap) {

        ArrayList<Double> pheno = new ArrayList<>();
        for (int i = 0; i < gtSamples.size(); i++) {
            Double valuetoAdd = phenoMap.getOrDefault(gtSamples.get(i), -9.0);
            if (!valuetoAdd.equals(-9.0)) {
                pheno.add(valuetoAdd);
            }
        }
        return mean(pheno.stream().mapToDouble(i -> i).toArray());
    }

    public static double percentile(ArrayList<Double> latencies, double percentile) {
        Collections.sort(latencies);
        double index = (double) Math.ceil(percentile / 100.0 * latencies.size());
        return latencies.get((int) (index - 1));
    }

}
