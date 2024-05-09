/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.mycompany.htssidra;

import java.util.ArrayList;
import java.util.Collections;

/**
 *
 * @author ealiyev
 */
public class Stat {

    private int ref_count_del;
    private int het_count_del;
    private int hom_count_del;
    private int ref_count_dup;
    private int het_count_dup;
    private int hom_count_dup;
    private int ref_count_inv;
    private int het_count_inv;
    private int hom_count_inv;
    private int gene_number_het_del;
    private int gene_number_het_dup;
    private int gene_number_het_inv;
    private int gene_number_hom_del;
    private int gene_number_hom_dup;
    private int gene_number_hom_inv;
    private int gene_omim_number_het_del;
    private int gene_omim_number_het_dup;
    private int gene_omim_number_het_inv;
    private int gene_omim_number_hom_del;
    private int gene_omim_number_hom_dup;
    private int gene_omim_number_hom_inv;
    private int pli_number_het_del;
    private int pli_number_het_dup;
    private int pli_number_het_inv;
    private int pli_number_hom_del;
    private int pli_number_hom_dup;
    private int pli_number_hom_inv;
    private ArrayList<Integer> del_het_sizes;
    private ArrayList<Integer> dup_het_sizes;
    private ArrayList<Integer> inv_het_sizes;
    private ArrayList<Integer> del_hom_sizes;
    private ArrayList<Integer> dup_hom_sizes;
    private ArrayList<Integer> inv_hom_sizes;

    public String calculateMean(ArrayList<Integer> sizes) {
        double sum = 0;
        for (int i = 0; i < sizes.size(); i++) {
            sum += sizes.get(i);
        }
        return String.valueOf(sum / sizes.size());
    }

    public String calculateMedian(ArrayList<Integer> sizes) {
        try {
            Collections.sort(sizes);
            int middle = sizes.size() / 2;
            if (sizes.size() % 2 == 1) {
                return String.valueOf(sizes.get(middle));
            } else {
                return String.valueOf((sizes.get(middle - 1) + sizes.get(middle)) / 2.0);
            }
        } catch (Exception e) {
            return "0";
        }

    }

    public String calculateTotalSize(ArrayList<Integer> sizes) {
        int sum = sizes.stream().mapToInt(Integer::intValue).sum();
        return String.valueOf(sum);
    }

    public String findMinimum(ArrayList<Integer> sizes) {
        if (sizes.size() > 0) {
            return String.valueOf(Collections.min(sizes));
        } else {
            return "NA";
        }
    }

    public String findMaximum(ArrayList<Integer> sizes) {
        if (sizes.size() > 0) {
            return String.valueOf(Collections.max(sizes));
        } else {
            return "NA";
        }
    }

    @Override
    public String toString() {
        return ref_count_del + "\t"
                + het_count_del + "\t"
                + hom_count_del + "\t"
                + ref_count_dup + "\t"
                + het_count_dup + "\t"
                + hom_count_dup + "\t"
                + ref_count_inv + "\t"
                + het_count_inv + "\t"
                + hom_count_inv + "\t"
                + gene_number_het_del + "\t"
                + gene_number_het_dup + "\t"
                + gene_number_het_inv + "\t"
                + gene_number_hom_del + "\t"
                + gene_number_hom_dup + "\t"
                + gene_number_hom_inv + "\t"
                + gene_omim_number_het_del + "\t"
                + gene_omim_number_het_dup + "\t"
                + gene_omim_number_het_inv + "\t"
                + gene_omim_number_hom_del + "\t"
                + gene_omim_number_hom_dup + "\t"
                + gene_omim_number_hom_inv + "\t"
                + pli_number_het_del + "\t"
                + pli_number_het_dup + "\t"
                + pli_number_het_inv + "\t"
                + pli_number_hom_del + "\t"
                + pli_number_hom_dup + "\t"
                + pli_number_hom_inv + "\t"
                + findMinimum(del_het_sizes) + "\t"
                + findMaximum(del_het_sizes) + "\t"
                + calculateMean(del_het_sizes) + "\t"
                + calculateMedian(del_het_sizes) + "\t"
                + calculateTotalSize(del_het_sizes) + "\t"
                + findMinimum(dup_het_sizes) + "\t"
                + findMaximum(dup_het_sizes) + "\t"
                + calculateMean(dup_het_sizes) + "\t"
                + calculateMedian(dup_het_sizes) + "\t"
                + calculateTotalSize(dup_het_sizes) + "\t"
                + findMinimum(inv_het_sizes) + "\t"
                + findMaximum(inv_het_sizes) + "\t"
                + calculateMean(inv_het_sizes) + "\t"
                + calculateMedian(inv_het_sizes) + "\t"
                + calculateTotalSize(inv_het_sizes) + "\t"
                + findMinimum(del_hom_sizes) + "\t"
                + findMaximum(del_hom_sizes) + "\t"
                + calculateMean(del_hom_sizes) + "\t"
                + calculateMedian(del_hom_sizes) + "\t"
                + calculateTotalSize(del_hom_sizes) + "\t"
                + findMinimum(dup_hom_sizes) + "\t"
                + findMaximum(dup_hom_sizes) + "\t"
                + calculateMean(dup_hom_sizes) + "\t"
                + calculateMedian(dup_hom_sizes) + "\t"
                + calculateTotalSize(dup_hom_sizes) + "\t"
                + findMinimum(inv_hom_sizes) + "\t"
                + findMaximum(inv_hom_sizes) + "\t"
                + calculateMean(inv_hom_sizes) + "\t"
                + calculateMedian(inv_hom_sizes) + "\t"
                + calculateTotalSize(inv_hom_sizes) + "\t";
    }

    public Stat() {
        ref_count_del = 0;
        het_count_del = 0;
        hom_count_del = 0;
        ref_count_dup = 0;
        het_count_dup = 0;
        hom_count_dup = 0;
        ref_count_inv = 0;
        het_count_inv = 0;
        hom_count_inv = 0;
        gene_number_het_del = 0;
        gene_number_het_dup = 0;
        gene_number_het_inv = 0;
        gene_number_hom_del = 0;
        gene_number_hom_dup = 0;
        gene_number_hom_inv = 0;
        gene_omim_number_het_del = 0;
        gene_omim_number_het_dup = 0;
        gene_omim_number_het_inv = 0;
        gene_omim_number_hom_del = 0;
        gene_omim_number_hom_dup = 0;
        gene_omim_number_hom_inv = 0;
        pli_number_het_del = 0;
        pli_number_het_dup = 0;
        pli_number_het_inv = 0;
        pli_number_hom_del = 0;
        pli_number_hom_dup = 0;
        pli_number_hom_inv = 0;
        del_het_sizes = new ArrayList<>();
        dup_het_sizes = new ArrayList<>();
        inv_het_sizes = new ArrayList<>();
        del_hom_sizes = new ArrayList<>();
        dup_hom_sizes = new ArrayList<>();
        inv_hom_sizes = new ArrayList<>();
    }

    /**
     * @return the ref_count_del
     */
    public int getRef_count_del() {
        return ref_count_del;
    }

    /**
     * @param ref_count_del the ref_count_del to set
     */
    public void setRef_count_del(int ref_count_del) {
        this.ref_count_del = ref_count_del;
    }

    /**
     * @return the het_count_del
     */
    public int getHet_count_del() {
        return het_count_del;
    }

    /**
     * @param het_count_del the het_count_del to set
     */
    public void setHet_count_del(int het_count_del) {
        this.het_count_del = het_count_del;
    }

    /**
     * @return the hom_count_del
     */
    public int getHom_count_del() {
        return hom_count_del;
    }

    /**
     * @param hom_count_del the hom_count_del to set
     */
    public void setHom_count_del(int hom_count_del) {
        this.hom_count_del = hom_count_del;
    }

    /**
     * @return the ref_count_dup
     */
    public int getRef_count_dup() {
        return ref_count_dup;
    }

    /**
     * @param ref_count_dup the ref_count_dup to set
     */
    public void setRef_count_dup(int ref_count_dup) {
        this.ref_count_dup = ref_count_dup;
    }

    /**
     * @return the het_count_dup
     */
    public int getHet_count_dup() {
        return het_count_dup;
    }

    /**
     * @param het_count_dup the het_count_dup to set
     */
    public void setHet_count_dup(int het_count_dup) {
        this.het_count_dup = het_count_dup;
    }

    /**
     * @return the hom_count_dup
     */
    public int getHom_count_dup() {
        return hom_count_dup;
    }

    /**
     * @param hom_count_dup the hom_count_dup to set
     */
    public void setHom_count_dup(int hom_count_dup) {
        this.hom_count_dup = hom_count_dup;
    }

    /**
     * @return the ref_count_inv
     */
    public int getRef_count_inv() {
        return ref_count_inv;
    }

    /**
     * @param ref_count_inv the ref_count_inv to set
     */
    public void setRef_count_inv(int ref_count_inv) {
        this.ref_count_inv = ref_count_inv;
    }

    /**
     * @return the het_count_inv
     */
    public int getHet_count_inv() {
        return het_count_inv;
    }

    /**
     * @param het_count_inv the het_count_inv to set
     */
    public void setHet_count_inv(int het_count_inv) {
        this.het_count_inv = het_count_inv;
    }

    /**
     * @return the hom_count_inv
     */
    public int getHom_count_inv() {
        return hom_count_inv;
    }

    /**
     * @param hom_count_inv the hom_count_inv to set
     */
    public void setHom_count_inv(int hom_count_inv) {
        this.hom_count_inv = hom_count_inv;
    }

    /**
     * @return the gene_number_het_del
     */
    public int getGene_number_het_del() {
        return gene_number_het_del;
    }

    /**
     * @param gene_number_het_del the gene_number_het_del to set
     */
    public void setGene_number_het_del(int gene_number_het_del) {
        this.gene_number_het_del = gene_number_het_del;
    }

    /**
     * @return the gene_number_het_dup
     */
    public int getGene_number_het_dup() {
        return gene_number_het_dup;
    }

    /**
     * @param gene_number_het_dup the gene_number_het_dup to set
     */
    public void setGene_number_het_dup(int gene_number_het_dup) {
        this.gene_number_het_dup = gene_number_het_dup;
    }

    /**
     * @return the gene_number_het_inv
     */
    public int getGene_number_het_inv() {
        return gene_number_het_inv;
    }

    /**
     * @param gene_number_het_inv the gene_number_het_inv to set
     */
    public void setGene_number_het_inv(int gene_number_het_inv) {
        this.gene_number_het_inv = gene_number_het_inv;
    }

    /**
     * @return the gene_number_hom_del
     */
    public int getGene_number_hom_del() {
        return gene_number_hom_del;
    }

    /**
     * @param gene_number_hom_del the gene_number_hom_del to set
     */
    public void setGene_number_hom_del(int gene_number_hom_del) {
        this.gene_number_hom_del = gene_number_hom_del;
    }

    /**
     * @return the gene_number_hom_dup
     */
    public int getGene_number_hom_dup() {
        return gene_number_hom_dup;
    }

    /**
     * @param gene_number_hom_dup the gene_number_hom_dup to set
     */
    public void setGene_number_hom_dup(int gene_number_hom_dup) {
        this.gene_number_hom_dup = gene_number_hom_dup;
    }

    /**
     * @return the gene_number_hom_inv
     */
    public int getGene_number_hom_inv() {
        return gene_number_hom_inv;
    }

    /**
     * @param gene_number_hom_inv the gene_number_hom_inv to set
     */
    public void setGene_number_hom_inv(int gene_number_hom_inv) {
        this.gene_number_hom_inv = gene_number_hom_inv;
    }

    /**
     * @return the gene_omim_number_het_del
     */
    public int getGene_omim_number_het_del() {
        return gene_omim_number_het_del;
    }

    /**
     * @param gene_omim_number_het_del the gene_omim_number_het_del to set
     */
    public void setGene_omim_number_het_del(int gene_omim_number_het_del) {
        this.gene_omim_number_het_del = gene_omim_number_het_del;
    }

    /**
     * @return the gene_omim_number_het_dup
     */
    public int getGene_omim_number_het_dup() {
        return gene_omim_number_het_dup;
    }

    /**
     * @param gene_omim_number_het_dup the gene_omim_number_het_dup to set
     */
    public void setGene_omim_number_het_dup(int gene_omim_number_het_dup) {
        this.gene_omim_number_het_dup = gene_omim_number_het_dup;
    }

    /**
     * @return the gene_omim_number_het_inv
     */
    public int getGene_omim_number_het_inv() {
        return gene_omim_number_het_inv;
    }

    /**
     * @param gene_omim_number_het_inv the gene_omim_number_het_inv to set
     */
    public void setGene_omim_number_het_inv(int gene_omim_number_het_inv) {
        this.gene_omim_number_het_inv = gene_omim_number_het_inv;
    }

    /**
     * @return the gene_omim_number_hom_del
     */
    public int getGene_omim_number_hom_del() {
        return gene_omim_number_hom_del;
    }

    /**
     * @param gene_omim_number_hom_del the gene_omim_number_hom_del to set
     */
    public void setGene_omim_number_hom_del(int gene_omim_number_hom_del) {
        this.gene_omim_number_hom_del = gene_omim_number_hom_del;
    }

    /**
     * @return the gene_omim_number_hom_dup
     */
    public int getGene_omim_number_hom_dup() {
        return gene_omim_number_hom_dup;
    }

    /**
     * @param gene_omim_number_hom_dup the gene_omim_number_hom_dup to set
     */
    public void setGene_omim_number_hom_dup(int gene_omim_number_hom_dup) {
        this.gene_omim_number_hom_dup = gene_omim_number_hom_dup;
    }

    /**
     * @return the gene_omim_number_hom_inv
     */
    public int getGene_omim_number_hom_inv() {
        return gene_omim_number_hom_inv;
    }

    /**
     * @param gene_omim_number_hom_inv the gene_omim_number_hom_inv to set
     */
    public void setGene_omim_number_hom_inv(int gene_omim_number_hom_inv) {
        this.gene_omim_number_hom_inv = gene_omim_number_hom_inv;
    }

    /**
     * @return the pli_number_het_del
     */
    public int getPli_number_het_del() {
        return pli_number_het_del;
    }

    /**
     * @param pli_number_het_del the pli_number_het_del to set
     */
    public void setPli_number_het_del(int pli_number_het_del) {
        this.pli_number_het_del = pli_number_het_del;
    }

    /**
     * @return the pli_number_het_dup
     */
    public int getPli_number_het_dup() {
        return pli_number_het_dup;
    }

    /**
     * @param pli_number_het_dup the pli_number_het_dup to set
     */
    public void setPli_number_het_dup(int pli_number_het_dup) {
        this.pli_number_het_dup = pli_number_het_dup;
    }

    /**
     * @return the pli_number_het_inv
     */
    public int getPli_number_het_inv() {
        return pli_number_het_inv;
    }

    /**
     * @param pli_number_het_inv the pli_number_het_inv to set
     */
    public void setPli_number_het_inv(int pli_number_het_inv) {
        this.pli_number_het_inv = pli_number_het_inv;
    }

    /**
     * @return the pli_number_hom_del
     */
    public int getPli_number_hom_del() {
        return pli_number_hom_del;
    }

    /**
     * @param pli_number_hom_del the pli_number_hom_del to set
     */
    public void setPli_number_hom_del(int pli_number_hom_del) {
        this.pli_number_hom_del = pli_number_hom_del;
    }

    /**
     * @return the pli_number_hom_dup
     */
    public int getPli_number_hom_dup() {
        return pli_number_hom_dup;
    }

    /**
     * @param pli_number_hom_dup the pli_number_hom_dup to set
     */
    public void setPli_number_hom_dup(int pli_number_hom_dup) {
        this.pli_number_hom_dup = pli_number_hom_dup;
    }

    /**
     * @return the pli_number_hom_inv
     */
    public int getPli_number_hom_inv() {
        return pli_number_hom_inv;
    }

    /**
     * @param pli_number_hom_inv the pli_number_hom_inv to set
     */
    public void setPli_number_hom_inv(int pli_number_hom_inv) {
        this.pli_number_hom_inv = pli_number_hom_inv;
    }

    /**
     * @return the del_het_sizes
     */
    public ArrayList<Integer> getDel_het_sizes() {
        return del_het_sizes;
    }

    /**
     * @param del_het_sizes the del_het_sizes to set
     */
    public void setDel_het_sizes(ArrayList<Integer> del_het_sizes) {
        this.del_het_sizes = del_het_sizes;
    }

    /**
     * @return the dup_het_sizes
     */
    public ArrayList<Integer> getDup_het_sizes() {
        return dup_het_sizes;
    }

    /**
     * @param dup_het_sizes the dup_het_sizes to set
     */
    public void setDup_het_sizes(ArrayList<Integer> dup_het_sizes) {
        this.dup_het_sizes = dup_het_sizes;
    }

    /**
     * @return the inv_het_sizes
     */
    public ArrayList<Integer> getInv_het_sizes() {
        return inv_het_sizes;
    }

    /**
     * @param inv_het_sizes the inv_het_sizes to set
     */
    public void setInv_het_sizes(ArrayList<Integer> inv_het_sizes) {
        this.inv_het_sizes = inv_het_sizes;
    }

    /**
     * @return the del_hom_sizes
     */
    public ArrayList<Integer> getDel_hom_sizes() {
        return del_hom_sizes;
    }

    /**
     * @param del_hom_sizes the del_hom_sizes to set
     */
    public void setDel_hom_sizes(ArrayList<Integer> del_hom_sizes) {
        this.del_hom_sizes = del_hom_sizes;
    }

    /**
     * @return the dup_hom_sizes
     */
    public ArrayList<Integer> getDup_hom_sizes() {
        return dup_hom_sizes;
    }

    /**
     * @param dup_hom_sizes the dup_hom_sizes to set
     */
    public void setDup_hom_sizes(ArrayList<Integer> dup_hom_sizes) {
        this.dup_hom_sizes = dup_hom_sizes;
    }

    /**
     * @return the inv_hom_sizes
     */
    public ArrayList<Integer> getInv_hom_sizes() {
        return inv_hom_sizes;
    }

    /**
     * @param inv_hom_sizes the inv_hom_sizes to set
     */
    public void setInv_hom_sizes(ArrayList<Integer> inv_hom_sizes) {
        this.inv_hom_sizes = inv_hom_sizes;
    }

}
