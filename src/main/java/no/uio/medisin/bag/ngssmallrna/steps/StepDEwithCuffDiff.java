/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import no.uio.medisin.bag.ngssmallrna.pipeline.ReferenceDataLocations;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * This step is not correctly written. There is no parsing of configuration data
 * to ensure values are correct. therefore there is no guarantee that it will 
 * run correctly.
 * 
 * @author joey
 */
public class StepDEwithCuffDiff extends NGSStep{
    
    static  Logger logger = LogManager.getLogger();

    public static final NGSStepSubclass STEP_SUBCLASS           = NGSStepSubclass.NOTDATABLE;

    public static final String STEP_ID_STRING = "RNASeqExpressionAnalysis";
    
    private static final String ID_PATH2TOPHAT = "pathToTopHat";
    private static final String ID_PATH2CUFFLINKS = "pathToCufflinks";
    private static final String ID_PATH2CUFFDIFF = "pathToCuffdiff";
    private static final String ID_PATH2CUFFMERGE = "pathToCuffmerge";
    
    private static final String ID_REF_GENOME = "referenceGenome";
    
    private static final String INFILE_EXTENSION_1 = "_1.fq.gz";
    private static final String INFILE_EXTENSION_2 = "_2.fq.gz";
    
    private static final String ID_THREADS = "noOfThreads";
    
    private String cufflinksSoftware = "";
    private String cuffdiffSoftware = "";
    private String cuffmergeSoftware = "";
    private String tophatSoftware = "";
    
    private String referenceGenomeIndex = "";
    
    private int noOfThreads = 1;
    
    
    public StepDEwithCuffDiff(){
        classSubtype = NGSStepSubclass.NOTDATABLE;
    }

    /**
     * 
     * @param sid 
     */
    public StepDEwithCuffDiff(InputDataForStep sid){
        classSubtype = NGSStepSubclass.NOTDATABLE;
        stepInputData = sid;
    }
    
    
    
    @Override
    public String shortStepDescription(){
      return "Perform differential expression analysis using CuffDiff.";
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "Perform differential expression analysis using CuffDiff.\n\n"
              + "NOTE: This step is not correctly written. "
              + "There is no parsing of configuration data to ensure values are "
              + "correct. therefore there is no guarantee that the step will run correctly";
        }
    

    @Override
    public void parseStepParameters() throws Exception{
      
    }
    


    /*
    There is no parsing of configuration data, so the step is prone to error
    */
    @Override
    public void parseConfigurationData(HashMap configData) throws Exception {
        logger.info(STEP_ID_STRING + ": verify configuation data");
        
        
    }


    @Override
    public void verifyInputData() throws IOException {
        logger.info(STEP_ID_STRING + ": verify input data");
        this.setPaths();
       
        
        //cheking tophat software
        if (new File(this.getTopHatSoftware()).exists() == false){
            logger.error("tophat software is not found at the location <" + this.getTopHatSoftware()+ "> " );
            throw new IOException("tophat software is not found at the location <" + this.getTopHatSoftware()+ "> ");
        }
        
        //cheking cufflinks software
        if (new File(this.getCufflinksSoftware()).exists() == false){
            logger.error("Cufflinks software is not found at the location <" + this.getCufflinksSoftware()+ "> " );
            throw new IOException("Cufflinks software is not found at the location <" + this.getCufflinksSoftware()+ "> ");
        }
        
        //cheking cuffdiff software
        if (new File(this.getCuffdiffSoftware()).exists() == false){
            logger.error("Cuffdiff software is not found at the location <" + this.getCuffdiffSoftware()+ "> " );
            throw new IOException("Cuffdiff software is not found at the location <" + this.getCuffdiffSoftware()+ "> ");
        }
        
        //cheking cuffmerge software
        if (new File(this.getCuffmergeSoftware()).exists() == false){
            logger.error("Cuffmerge software is not found at the location <" + this.getCuffmergeSoftware()+ "> " );
            throw new IOException("Cuffmerge software is not found at the location <" + this.getCuffmergeSoftware()+ "> ");
        }
        
      
        
        //checking bowtie genome index
        String pathToBowtieGenomeIndex = this.cleanPath(getStepInputData().getDataLocations().getGenomeRootFolder()
                + FILESEPARATOR +this.getReferenceGenomeIndex() + FILESEPARATOR + ReferenceDataLocations.ID_REL_BOWTIE_PATH + ".1.ebwt");
        if(new File(pathToBowtieGenomeIndex).exists() == false){
            logger.error("bowtie genome index <" + pathToBowtieGenomeIndex + "> not found");
            throw new IOException("bowtie genome index <" + pathToBowtieGenomeIndex + "> not found");
        }
        
        //checking for the input data file
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        while(itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            String fastaFile1 = (String)sampleData.getFastqFile1();
            String fastaFile2 = (String)sampleData.getFastqFile2();
//            String fastaFile3 = (String)sampleData.getFastqFile3();
            
            //checking for fastFile1
            if (sampleData.getFastqFile1() == null){
                logger.error("no forward fastq file is specified");
                throw new IOException("no forward fastq file is specified");
            }
            else{
                if(fastaFile1.endsWith(".fastq" ) || fastaFile1.endsWith(".fq") || fastaFile1.endsWith(".fasta")){
                    if (fastaFile1.contains("_1")){
                        if(new File(this.cleanPath(fastaFile1)).exists() == false){
                            logger.error("The forward reads fastFile1<" + fastaFile1 + "> does not exist");
                            throw new IOException("The forward reads fastFile1<" + fastaFile1 + "> does not exist");
                        }
                    }
                    else{
                        logger.error("this fastaFile1 must be specified as forward sequence with _1");
                        throw new IOException("this fastaFile1 must be specified as forward sequence with _1");
                    } 
                }
                else{                    
                    logger.error("this data set must be in fastq or fasta format");
                    throw new IOException("this dataset must be in fastq or fasta format");
                }
            }
            
            //checking for fastFile2
            if (sampleData.getFastqFile2() == null){
                logger.error("no reverse fastq file is specified");
                throw new IOException("no reverse fastq file is specified");
            }
            else{
                if(fastaFile2.endsWith(".fastq" ) || fastaFile2.endsWith(".fq") || fastaFile2.endsWith(".fasta")){
                    if (fastaFile2.contains("_2")){
                        if(new File(this.cleanPath(fastaFile2)).exists() == false){
                            logger.error("The reverse reads fastaFile2<" + fastaFile2 + "> does not exist");
                            throw new IOException("The reverse reads fastaFile2<" + fastaFile2 + "> does not exist");
                        }
                    }
                    else{
                        logger.error("this fastaFile2 must be specified as reverse sequence with _2");
                        throw new IOException("this fastaFile2 must be specified as reverse sequence with _2");
                    } 
                }
                else{
                    
                    logger.error("this data set must be in fastq or fasta format");
                    throw new IOException("this dataset must be in fastq or fasta format");
                }
            }
            
//            //checking for fastFile3 ---- additional unpaired reads
//            if (sampleData.getFastqFile3() == null){
//                logger.info("no additional unpaired reads is specified");
//            }
//            else{
//                if(fastaFile3.endsWith(".fastq" ) || fastaFile3.endsWith(".fq") || fastaFile3.endsWith(".fasta")){
//                    if (fastaFile3.contains("_1") || fastaFile3.contains("_2")){
//                        if(new File(this.cleanPath(fastaFile1)).exists() == false){
//                            logger.error("The additional unpaired reads fastaFile3<" + fastaFile3 + "> does not exist");
//                            throw new IOException("The additional unpaired reads fastaFile3<" + fastaFile3 + "> does not exist");
//                        }
//                    }
//                    else{
//                        logger.error("this fastaFile3 is additional unpaired reads, should not contain _1 or _2");
//                        throw new IOException("this fastaFile3 is additional unpaired reads, should not contain _1 or _2");
//                    } 
//                }
//                else{                   
//                    logger.error("this data set must be in fastq or fasta format");
//                    throw new IOException("this dataset must be in fastq or fasta format");
//                }
//            }
            
            logger.info("passed");
        }  
    }

    @Override
    public void verifyOutputData()  {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public NGSStepSubclass getStepSubclass(){
        return STEP_SUBCLASS;
    }
    
    
    @Override
    public HashMap generateExampleConfigurationData() {
        logger.info(STEP_ID_STRING + ": generate example configuration data");

        HashMap configData = new HashMap();
        HashMap paramData = new HashMap();

        paramData.put(ID_PATH2TOPHAT, "/usr/local/bin/tophat");
        paramData.put(ID_PATH2CUFFLINKS, "/usr/local/bin/cufflinks");
        paramData.put(ID_PATH2CUFFDIFF, "/usr/local/bin/cuffdiff");
        paramData.put(ID_PATH2CUFFMERGE, "/usr/local/bin/cuffmerge");
        paramData.put(ID_REF_GENOME, "/data/ngsdata/genome");
        paramData.put(ID_THREADS, 4);
        
        configData.put(STEP_ID_STRING, paramData);

        return configData;
    }

    
    
    
    @Override
    public void execute() throws IOException {
        logger.info(STEP_ID_STRING + ":exexute");
        Boolean fBoolean = new File(outFolder).mkdir();
        if(fBoolean) {
            logger.info("created output folder <" + outFolder + "> for results");
        }
        
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        while(itSD.hasNext()){
            ArrayList<String> cmd = new ArrayList<String>();
            //todo
//            try {
//                SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
//                
//            } catch (IOException | InterruptedException ex) {
//                logger.error("error executing expression  command\n");
//                logger.error(cmd);
//                logger.error(ex.toString());
//                throw new IOException(STEP_ID_STRING + ": \"error executing Bowtie Mapping command " + cmd);
//            }

        }
    }

    /**
     * 
     * @return tophatSoftware
     */
    private String getTopHatSoftware() {
        return this.tophatSoftware;
    }
    
    /**
     * set
     * @param tophatsoftware 
     */
    private void setTopHatSoftware(String tophatsoftware){
        this.tophatSoftware = tophatsoftware;
    }

    /**
     * 
     * @return cufflinksSoftware
     */
    private String getCufflinksSoftware() {
        return this.cufflinksSoftware;
    }
    
    /**
     * set
     * @param cufflinkssoftware 
     */
    private void setCufflinksSoftware(String cufflinkssoftware){
        this.cufflinksSoftware = cufflinkssoftware;
    }

    /**
     * 
     * @return cuffdiffSoftware
     */
    private String getCuffdiffSoftware() {
        return this.cuffdiffSoftware;
    }
    
    /**
     * set
     * @param cuffdiffsoftware 
     */
    private void setCuffdiffSoftware(String cuffdiffsoftware){
        this.cuffdiffSoftware = cuffdiffsoftware;
    }

    /**
     * 
     * @return referenceGenomeIndex
     */
    private String getReferenceGenomeIndex() {
        return this.referenceGenomeIndex;
    }
    
    /**
     * set
     * @param referenceGenomeIndexString 
     */
    private void setReferenceGenomeIndex(String referenceGenomeIndexString){
        this.referenceGenomeIndex = referenceGenomeIndexString;
    }

    /**
     * 
     * @return cuffmergeSoftware
     */
    private String getCuffmergeSoftware() {
        return this.cuffmergeSoftware;
    }
    
    /**
     * set
     * @param cuffmergesoftware 
     */
    private void setCuffmergeSoftware(String cuffmergesoftware){
        this.cuffmergeSoftware = cuffmergesoftware;
    }
   
}
