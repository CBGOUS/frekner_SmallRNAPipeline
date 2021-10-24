/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

import org.apache.commons.lang3.Validate;

/**
 * stores information for a sample entry specified in a run data file
 * handles both single and paired end reads
 * for single end reads, the second file will be set to null
 * the entries correspond to the columns in a tab delimited data file that
 * is specified as part of a run. While fastqFile1 and fastqFile2 correspond 
 * to the filenamesfor the two fastq files, the type of information contained
 * in the other variables can be a little varied:
 *   
 *  - datasource: typically something such as a host name
 *  - condition : generally used for differential expression analysis 
 *                (e.g. healthy vs disease, different cell type...)
 *  - time      : for time course experiments
 *  - notes     : any general information about the sample which will not be
 *                used in the analysis
 * 
 * @author sr
 */
public class SampleDataEntry {
    
    private String fastqFile1;
    private String fastqFile2;
    private String dataSource;
    private String condition;
    private String time;
    private String note;

    
    
    public SampleDataEntry(String fname1, String fname2, String source, String cond, String t, String n){
        
        fastqFile1 = fname1;
        fastqFile2 = fname2;
        dataSource = source;
        condition = cond;
        time = t;
        note = n;
        
        Validate.notNull(fastqFile1, "fastq file 1 must be specified");
        Validate.notNull(dataSource, "source must be specified");
        
    }
    
    
    /*
        @return the instance as a string
    */
    @Override
    public String toString(){
        String s = "fastqFile1:\t" + fastqFile1 + "\t"
                + "fastqFile2:\t" + fastqFile2 + "\t"
                + "dataSource:\t" + dataSource + "\t"
                + "condition:\t" + condition + "\t"
                + "time:\t" + time + "\t"
                + "note:\t" + note;
                
        return s;
    };
    
    
    /**
     * @return the dataFile
     */
    public String getFastqFile1() {
        return fastqFile1;
    }
    

    /**
     * @param dataFile the dataFile to set
     */
    public void setDataFile(String dataFile) {
        this.fastqFile1 = dataFile;
    }

    /**
     * @return the dataSource
     */
    public String getDataSource() {
        return dataSource;
    }

    /**
     * @param dataSource the dataSource to set
     */
    public void setDataSource(String dataSource) {
        this.dataSource = dataSource;
    }

    /**
     * @return the Condition
     */
    public String getCondition() {
        return condition;
    }

    /**
     * @param Condition the Condition to set
     */
    public void setCondition(String Condition) {
        this.condition = Condition;
    }

    /**
     * @return the time
     */
    public String getTime() {
        return time;
    }

    /**
     * @param time the time to set
     */
    public void setTime(String time) {
        this.time = time;
    }

    /**
     * @return the note
     */
    public String getNote() {
        return note;
    }

    /**
     * @param note the note to set
     */
    public void setNote(String note) {
        this.note = note;
    }

    /**
     * @return the fastqFile2
     */
    public String getFastqFile2() {
        return fastqFile2;
    }
    

    /**
     * @param fastqFile2 the fastqFile2 to set
     */
    public void setFastqFile2(String fastqFile2) {
        this.fastqFile2 = fastqFile2;
    }
    
}
