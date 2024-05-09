/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.mycompany.htssidra;

import java.io.IOException;

/**
 *
 * @author ealiyev
 */
public class Main {

    public static void main(String[] args) throws IOException, InterruptedException {

        if (args.length < 1) {
            System.out.println("Wrong input");
        } else {
            switch (args[0]) {
                case "mergephenos":
                    Helper.mergePhenotypes(args[1], args[2]);
                    break;
                case "convert":
                    Helper.calculateSVTypeSvaba(args[1]);
                    break;
                case "merge":
                    Helper.mergeFiles("C:\\Users\\ealiyev\\Dropbox\\SV Project\\AnnotatedSvabaPMC\\Homo");
                    break;
                case "cnvtopredictor":
                    Helper.popVCFtoPoppante(args[1]);
                    break;
                case "renamesamplename":
                    Helper.renameSampleName(args[1]);
                    break;
                case "extractnamemapping":
                    Helper.extractSampleBam(args[1]);
                    break;
                case "processAnnotatedFile":
                    Helper.processAnnotatedFile(args[1]);
                    break;
                case "processAnnotatedFileAnnotSV_3_8":
                    Helper.processAnnotatedFile_3_8(args[1]);
                    break;
                case "processAnnotatedFileSplit":
                    Helper.processAnnotatedFileSplit(args[1]);
                    break;
                case "countgtmultisampleVCF":
                    Helper.vcfStats(args[1], args[2], args[3]);
                    break;
                case "createDenovoHomozygousReport":
                    Helper.generateReportsbyPED(args[1], args[2], args[3]);
                    break;
                case "annotateControlCaseCounts":
                    Helper.phenotypeAnnotation(args[1], args[2]);
                    break;
                case "annotatewithSubpops":
                    Helper.annotateFilewithSubPopFrequencies(args[1], args[2]);
                    break;
                case "ranker":
                    Helper.ranker(args[1], args[2]);
                    break;
                case "dummyVCF":
                    Helper.dummyVCF(args[1], args[2]);
                    break;
                case "svInspector":
                    Helper.svInspector(args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8]);
                    break;
                case "svClusterMCNV":
                    Helper.clusterVCF(args[1]);
                    break;
                case "annotateVCFwithPhenotypes":
                    Helper.annotateVCFwithPhenos(args[1], args[2]);
                    break;
                case "adjustGenotype":
                    Helper.adjustGenotype(args[1]);
                    break;
                case "adjustALT":
                    Helper.adjustALT(args[1]);
                    break;

            }
        }

    }
}
