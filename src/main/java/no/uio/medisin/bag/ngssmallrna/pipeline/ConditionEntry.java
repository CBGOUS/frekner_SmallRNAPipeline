/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

/**
 * This stores a condition entry as defined in the TSV data file specified 
 * at run time. e.g. if the TSV file has the format
 *     File                Source    Condition  Time	Note
 *   2726_S3.fastq          Msy        ph5.5     t0      Msy
 *   2726_S4.fastq          Msy        ph6.1     t0      Msy
 *
 * then there are two conditions, "ph5.5" & "ph6.1"
 * 
 * this can is used to consolidate data from different samples but identical conditions
 * 
 * @author sr
 */
public class ConditionEntry {
    private String name;
    private String dataset1;
    private String dataset2;
    
    public ConditionEntry(String n, String d1, String d2){
        
        name = n;
        dataset1 = d1;
        dataset2 = d2;
        
    }

    /**
     * print this instance
     * @return 
     */
    public String printEntry(){
        return name + "\t" + dataset1 + "\t" + dataset2;
    }
    
    
    /**
     * @return the name
     */
    public String getName() {
        return name;
    }

    /**
     * @return the dataset1
     */
    public String getDataset1() {
        return dataset1;
    }

    /**
     * @return the dataset2
     */
    public String getDataset2() {
        return dataset2;
    }

    /**
     * @param dataset1 the dataset1 to set
     */
    public void setDataset1(String dataset1) {
        this.dataset1 = dataset1;
    }

    /**
     * @param dataset2 the dataset2 to set
     */
    public void setDataset2(String dataset2) {
        this.dataset2 = dataset2;
    }
    

    
}
