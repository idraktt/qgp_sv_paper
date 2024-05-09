/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.mycompany.htssidra;

import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author ealiyev
 */
public class Phenotype {

    private String name;
    private HashMap<String, Double> pheno_raw;
    private HashMap<String, Double> pheno_filtered;
    private HashMap<String, Double> pheno_Rank;

    public Phenotype(String name, HashMap<String, Double> pheno_raw, HashMap<String, Double> pheno_filtered, HashMap<String, Double> pheno_Rank) {
        this.name = name;
        this.pheno_raw = pheno_raw;
        this.pheno_filtered = pheno_filtered;
        this.pheno_Rank = pheno_Rank;
    }

    /**
     * @return the name
     */
    public String getName() {
        return name;
    }

    /**
     * @param name the name to set
     */
    public void setName(String name) {
        this.name = name;
    }

    /**
     * @return the pheno_raw
     */
    public HashMap<String, Double> getPheno_raw() {
        return pheno_raw;
    }

    /**
     * @param pheno_raw the pheno_raw to set
     */
    public void setPheno_raw(HashMap<String, Double> pheno_raw) {
        this.pheno_raw = pheno_raw;
    }

    /**
     * @return the pheno_filtered
     */
    public HashMap<String, Double> getPheno_filtered() {
        return pheno_filtered;
    }

    /**
     * @param pheno_filtered the pheno_filtered to set
     */
    public void setPheno_filtered(HashMap<String, Double> pheno_filtered) {
        this.pheno_filtered = pheno_filtered;
    }

    /**
     * @return the pheno_Rank
     */
    public HashMap<String, Double> getPheno_Rank() {
        return pheno_Rank;
    }

    /**
     * @param pheno_Rank the pheno_Rank to set
     */
    public void setPheno_Rank(HashMap<String, Double> pheno_Rank) {
        this.pheno_Rank = pheno_Rank;
    }

}
