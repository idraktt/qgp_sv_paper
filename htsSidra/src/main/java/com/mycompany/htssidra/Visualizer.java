/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.mycompany.htssidra;

import au.com.bytecode.opencsv.CSVReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.lang3.StringUtils;

/**
 *
 * @author ealiyev
 */
public class Visualizer {

  

    public String getBamFileNamebySampleName() {
        return "";
    }

    public void compareVectors(String firstVector, String secondVector) {
        StringUtils.difference(firstVector, firstVector);

        int numberofMatches = StringUtils.countMatches(firstVector, firstVector);

    }

}
