/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import no.uio.medisin.bag.ngssmallrna.pipeline.ReferenceDataLocations;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Perform TopHat mapping of RNASeq reads
 * @author joey
 */
public class StepTopHatAlignRNAReads extends NGSStep{
    
    static  Logger logger = LogManager.getLogger();
    
    public static final NGSStepSubclass STEP_SUBCLASS           = NGSStepSubclass.DATABLE;    

    public static final String STEP_ID_STRING = "TopHatAlignRNAReads";
    
    private static final String ID_PATH2BOWTIE = "pathToBowtie";
    private static final String ID_APTH2TOPHAT = "pathToTopHat";
    
    private static final String ID_REF_GENOME = "referenceGenome";
   
    private static final String ID_READS_MODEL = "readsModel";
    private static final String ID_MIXED_READS = "mixedReads";
    
    private static final String ID_READ_MISMATCHES = "readMismatches";
    private static final String ID_READ_GAP_LEN = "readGapLength";
    private static final String ID_READ_EDIT_DIST = "readEditDist";
    private static final String ID_READ_REALIGN_EDIT_DIST = "readRealignEditDist";
    
    private static final String INFILE_EXTENSION_1 = "_1.fq.gz";
    private static final String INFILE_EXTENSION_2 = "_2.fq.gz";
    
    private static final String ID_USE_BOWTIE1 = "useBowtie1";
    
    private static final String ID_THREADS = "noOfThreads";
    
    private String bowtieSoftware = "";
    private String tophatSoftware = "";
    private String referenceGenomeIndex = "";
    private String readModel = "";
    private String mixedReads = "";
    private int readMismatches = 2; 
    private int readGapLength = 2;
    private int readEditDist = 2;
    private int readRealignEditDist = 2;      
    private String useBowtie1 = "";       
    private int noOfThreads = 1;
//    private boolean additionalSEReads = false;
    
    
    
    public StepTopHatAlignRNAReads(){
        classSubtype = NGSStepSubclass.DATABLE;
    }

    /**
     * 
     * @param sid stepInputData
     */
    public StepTopHatAlignRNAReads(InputDataForStep sid){
        classSubtype = NGSStepSubclass.DATABLE;
        stepInputData = sid;
    }


    @Override
    public String shortStepDescription(){
      return "Perform TopHat mapping of RNASeq read";
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "Perform TopHat mapping of RNASeq read.\n"
              + "The step requires the TopHat software to be installed.\n"
              + "The location can be specified in the YAML configuration file\n"
              + "together with additional parameters\n";
    }


    @Override
    public void parseStepParameters() throws Exception{
      
    }
    
    
    
    @Override
    public void verifyInputData() throws IOException {
        logger.info(STEP_ID_STRING + ": verify input data");
        this.setPaths();
        
        //checking bowtie software
        if(new File(this.getBowtieSoftware()).exists() == false){
            logger.error("bowtie software is not found at the location <" + this.getBowtieSoftware()+ "> " );
            throw new IOException("bowtie software is not found at the location <" + this.getBowtieSoftware()+ "> " );
        }
        
        //cheking tophat software
        if (new File(this.getTopHatSoftware()).exists() == false){
            logger.error("tophat software is not found at the location <" + this.getTopHatSoftware()+ "> " );
            throw new IOException("tophat software is not found at the location <" + this.getTopHatSoftware()+ "> ");
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
            
            //checking for fastFile3 ---- additional unpaired reads
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
    
    
    
    
    
    /**
     * simply checking configuration file has all entries we need 
     * @param configData
     * @throws Exception 
     */
    @Override
    public void parseConfigurationData(HashMap configData) throws Exception {
        logger.info(STEP_ID_STRING + ": verify configuation data");
        
        String[] id_itmes = {ID_PATH2BOWTIE, ID_APTH2TOPHAT, ID_REF_GENOME, ID_READS_MODEL, ID_MIXED_READS,
        ID_READ_MISMATCHES, ID_READ_GAP_LEN, ID_READ_EDIT_DIST, ID_READ_REALIGN_EDIT_DIST,
        INFILE_EXTENSION_1, INFILE_EXTENSION_2, ID_USE_BOWTIE1, ID_THREADS};
        
        for (int i =0; i <= id_itmes.length-1; i++){
            if(configData.get(id_itmes[i]) == null){
                logger.error("<" + configData.get(id_itmes[i]) + "> : Missing definition in configutation file");
                throw new NullPointerException("<" + configData.get(id_itmes[i]) + "> : Missing definition in configutation file");
            }
        }
        
        //checking NoOfThreads
        try{
            this.setNoOfThreads((Integer) configData.get(ID_THREADS));
        }
        catch(NumberFormatException exNm){
            logger.error(ID_THREADS + " <" + ID_THREADS + "> is not an integer");
            throw new NumberFormatException(ID_THREADS + " <" + ID_THREADS + "> is not an integer");
        }        
        if (this.getNoOfThreads() <= 0){
            logger.error(ID_THREADS + " <" + configData.get(ID_THREADS) + "> must be positive integer");
            throw new IllegalArgumentException(ID_THREADS + " <" + configData.get(ID_THREADS) + "> must be positive integer");
        }
        
        //checking ID_READ_GAP_LEN
        try{
            this.setReadGapLength((Integer) configData.get(ID_READ_GAP_LEN));
        }
        catch(NumberFormatException exNm){
            logger.error(ID_READ_GAP_LEN + " <" + ID_READ_GAP_LEN + "> is not an integer");
            throw new NumberFormatException(ID_READ_GAP_LEN + " <" + ID_READ_GAP_LEN + "> is not an integer");
        }        
        if (this.getReadGapLength() <= 0){
            logger.error(ID_READ_GAP_LEN + " <" + configData.get(ID_READ_GAP_LEN) + "> must be positive integer");
            throw new IllegalArgumentException(ID_READ_GAP_LEN + " <" + configData.get(ID_READ_GAP_LEN) + "> must be positive integer");
        }
        
        //checking ID_READ_MISMATCHES
        try{
            this.setNoOfReadMismatches((Integer) configData.get(ID_READ_MISMATCHES));
        }
        catch(NumberFormatException exNm){
            logger.error(ID_READ_MISMATCHES + " <" + ID_READ_MISMATCHES + "> is not an integer");
            throw new NumberFormatException(ID_READ_MISMATCHES + " <" + ID_READ_MISMATCHES + "> is not an integer");
        }        
        if (this.getNoOfReadMismatches() <= 0){
            logger.error(ID_READ_MISMATCHES + " <" + configData.get(ID_READ_MISMATCHES) + "> must be positive integer");
            throw new IllegalArgumentException(ID_READ_MISMATCHES + " <" + configData.get(ID_READ_MISMATCHES) + "> must be positive integer");
        }
        
        //checking ID_READ_EDIT_DIST
        try{
            this.setReadEditDist((Integer) configData.get(ID_READ_EDIT_DIST));
        }
        catch(NumberFormatException exNm){
            logger.error(ID_READ_EDIT_DIST + " <" + ID_READ_EDIT_DIST + "> is not an integer");
            throw new NumberFormatException(ID_READ_EDIT_DIST + " <" + ID_READ_EDIT_DIST + "> is not an integer");
        }        
        if (this.getReadEditDist() <= 0){
            logger.error(ID_READ_EDIT_DIST + " <" + configData.get(ID_READ_EDIT_DIST) + "> must be positive integer");
            throw new IllegalArgumentException(ID_READ_EDIT_DIST + " <" + configData.get(ID_READ_EDIT_DIST) + "> must be positive integer");
        }
        
        //checking ID_READ_REALIGN_EDIT_DIST
        try{
            this.setReadRealignEditDist((Integer) configData.get(ID_READ_REALIGN_EDIT_DIST));
        }
        catch(NumberFormatException exNm){
            logger.error(ID_READ_REALIGN_EDIT_DIST + " <" + ID_READ_REALIGN_EDIT_DIST + "> is not an integer");
            throw new NumberFormatException(ID_READ_REALIGN_EDIT_DIST + " <" + ID_READ_REALIGN_EDIT_DIST + "> is not an integer");
        }        
        if (this.getReadRealignEditDist() <= 0){
            logger.error(ID_READ_REALIGN_EDIT_DIST + " <" + configData.get(ID_READ_REALIGN_EDIT_DIST) + "> must be positive integer");
            throw new IllegalArgumentException(ID_READ_REALIGN_EDIT_DIST + " <" + configData.get(ID_READ_REALIGN_EDIT_DIST) + "> must be positive integer");
        }
        
        //checking ID_READS_MODEL
        if(configData.get(ID_READS_MODEL).toString().contains("SE") || configData.get(ID_READS_MODEL).toString().contains("SE")){
            this.setReadModel((String)configData.get(ID_READS_MODEL));
        }
        else{
            logger.error(ID_READS_MODEL + "<" +ID_READS_MODEL + "> must be set as \"SE\" or \"PE\". \"SE\" for single end, \"PE\" fpr paried end");
            throw new IllegalArgumentException(ID_READS_MODEL + "<" +ID_READS_MODEL + "> must be set as \"SE\" or \"PE\". \"SE\" for single end, \"PE\" fpr paried end");
        }
                
        
        //checking ID_MIXED_READS
        if(configData.get(ID_MIXED_READS).toString().contains("true") || configData.get(ID_MIXED_READS).toString().contains("false")){
            this.setMixedReads((String)configData.get(ID_MIXED_READS));
        }
        else{
            logger.error(ID_READS_MODEL + "<" +ID_READS_MODEL + "> must be set as \"true\" or \"false\"");
            throw new IllegalArgumentException(ID_READS_MODEL + "<" +ID_READS_MODEL + "> must be set as \"true\" or \"false\"");
        }
        
        
        this.setReferenceGenomeIndex((String)configData.get(ID_REF_GENOME));
        if(this.getReferenceGenomeIndex().length() != 3){
            logger.error(ID_REF_GENOME + "<" + configData.get(ID_REF_GENOME) +"> must be a 3 letter string");
            throw new IllegalArgumentException(ID_REF_GENOME + "<" + configData.get(ID_REF_GENOME) +"> must be a 3 letter string");
        }
        
        logger.info("passed");
        
  

    }

    @Override
    public HashMap generateExampleConfigurationData() {
        logger.info(STEP_ID_STRING + ": generate example configuration data");

        HashMap configData = new HashMap();
        HashMap paramData = new HashMap();

        paramData.put(ID_PATH2BOWTIE, "/usr/local/bin/bowtie");
        paramData.put(ID_APTH2TOPHAT, "/usr/local/bin/tophat");
        paramData.put(ID_REF_GENOME, "/data/ngsData/genome/hsa");
        paramData.put(ID_READS_MODEL, "SE");
        paramData.put(ID_MIXED_READS, "true");
        paramData.put(ID_READ_MISMATCHES, 2);
        paramData.put(ID_READ_GAP_LEN, 2);
        paramData.put(ID_READ_EDIT_DIST, 2);
        paramData.put(ID_READ_REALIGN_EDIT_DIST, 2);
        paramData.put(ID_USE_BOWTIE1, "true");
        paramData.put(ID_THREADS, 4);
        
        configData.put(STEP_ID_STRING, paramData);

        return configData;
    }

    /**
     * run the tophat command
     * @throws IOException 
     */
    @Override
    public void execute() throws IOException {
        this.setPaths();
        this.verifyInputData();
        
        Boolean outFolderBoolean = new File(outFolder).mkdir();
        if (outFolderBoolean){
            logger.info("created output folder <" + outFolder + "> for results");
        }
        
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        String cmdTophat = "";
        while (itSD.hasNext()){
            try{
                SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
                ArrayList<String> tophatCMD = new ArrayList<>();
                
                tophatCMD.add("tophat");
                //-N/--read-mismatches	Final read alignments having more than these many mismatches are discarded. The default is 2.
                tophatCMD.add("-N " + this.getNoOfReadMismatches());
                //Final read alignments having more than these many total length of gaps are discarded. The default is 2.
                tophatCMD.add("--read-gap-length " + this.getReadGapLength());
                //Final read alignments having more than these many edit distance are discarded. The default is 2.
                tophatCMD.add("--read-edit-dist" + this.getReadEditDist());
                //The default value for this option is set such that TopHat will not try to realign reads already mapped in earlier steps.
                tophatCMD.add("--read-realign-edit-dist " + this.getReadRealignEditDist());
                tophatCMD.add("-p "+ this.getNoOfThreads());
                //Uses Bowtie1 instead of Bowtie2. If you use colorspace reads, 
                //you need to use this option as Bowtie2 does not support colorspace reads.
                if(useBowtie1.equals("true")){
                    tophatCMD.add("--bowtie1");
                }
                
                tophatCMD.add(referenceGenomeIndex);
                tophatCMD.add(inFolder + FILESEPARATOR + sampleData.getFastqFile1());
                tophatCMD.add(inFolder + FILESEPARATOR + sampleData.getFastqFile2());
                
//                if(sampleData.getFastqFile3() != null){
//                    tophatCMD.add(inFolder + FILESEPARATOR + sampleData.getFastqFile3());
//                }
                
                cmdTophat = this.cleanPath(StringUtils.join(tophatCMD, " "));
                logger.info("TopHat run command:\t" + cmdTophat + "\nProcessing input files:"+ sampleData.getFastqFile1()+ " and " + sampleData.getFastqFile2());

                Runtime rt = Runtime.getRuntime();
                Process proc = rt.exec(cmdTophat);
                BufferedReader brStdin  = new BufferedReader(new InputStreamReader(proc.getInputStream()));
                BufferedReader brStdErr = new BufferedReader(new InputStreamReader(proc.getErrorStream()));
                
                String line = null;
                    logger.info("<OUTPUT>");
                    while ( (line = brStdin.readLine()) != null)
                        logger.info(line);
                    logger.info("</OUTPUT>");

                    logger.info("<ERROR>");
                    while ( (line = brStdErr.readLine()) != null)
                        logger.info(line);
                    logger.info("</ERROR>");
                
                
                int exitVal = proc.waitFor();            
                logger.info("Process exitValue: " + exitVal);   
                
                brStdin.close();
                brStdErr.close();
            }
            catch(IOException| InterruptedException ex) {
                logger.error("error executing TopHat command\n" + ex.toString());
                throw new IOException(STEP_ID_STRING + "error executing TopHat command " + cmdTophat);
            } 
        }
        logger.info(STEP_ID_STRING + ": completed");
        
            
        /**
         * tophat [options]* <genome_index_base> PE_reads_1.fq.gz,SE_reads.fa PE_reads_2.fq.gz
         */
        
        
                
    }

    /**
     * 
     * @return bowtieSoftware
     */    
    private String getBowtieSoftware() {
        return bowtieSoftware;
    }
    
    /**
     * 
     * @param bowtieSoftware set
     */
    private void setBowtieSoftware(String bowtieSoftware) {
        this.bowtieSoftware = bowtieSoftware;
    }
    
    /**
     * 
     * @return tophatSoftware
     */
    private String getTopHatSoftware() {
        return tophatSoftware;
    }
    
    /**
     * 
     * @param tophatSoftwatr set
     */
    private void setTopHatSoftware(String tophatSoftwatr) {
        this.tophatSoftware = tophatSoftwatr;
    }

    /**
     * 
     * @return referenceGenomeIndex
     */
    private String getReferenceGenomeIndex() {
        return referenceGenomeIndex;
    }
    
    /**
     * 
     * @param referenceGenomeIndex  set
     */
    private void setReferenceGenomeIndex(String referenceGenomeIndex){
        this.referenceGenomeIndex = referenceGenomeIndex;
    }

    /**
     * 
     * @return noOfThreads
     */
    private int getNoOfThreads() {
        return  noOfThreads;
    }
    
    /**
     * 
     * @param intThreads set
     */
    private void setNoOfThreads(Integer intThreads) {
        this.noOfThreads = intThreads;
    }

    /**
     * 
     * @param readMismatches set
     */
    private void setNoOfReadMismatches(Integer intReadMismatches) {
        this.readMismatches = intReadMismatches;
    }

    /**
     * 
     * @return readMismatches
     */
    private int getNoOfReadMismatches() {
        return readMismatches;
    }

    /**
     * 
     * @return readGapLength
     */
    private int getReadGapLength() {
        return readGapLength;
    }

    /**
     * 
     * @param intReadGapLength set
     */
    private void setReadGapLength(Integer intReadGapLength) {
        this.readGapLength = intReadGapLength;
    }

    /**
     * 
     * @return readEditDist
     */
    private int getReadEditDist() {
        return readEditDist;
    }

    /**
     * 
     * @param intReadEditDist set
     */
    private void setReadEditDist(Integer intReadEditDist) {
        this.readEditDist = intReadEditDist;
    }

    /**
     * 
     * @return readRealignEditDist
     */
    private int getReadRealignEditDist() {
        return readRealignEditDist;
    }

    /**
     * 
     * @param intReadRealignEditDist set 
     */
    private void setReadRealignEditDist(Integer intReadRealignEditDist) {
        this.readRealignEditDist = intReadRealignEditDist;
    }
    
   

    /**
     * 
     * @param strReadModel set
     */
    private void setReadModel(String strReadModel) {
        this.readModel = strReadModel;
    }
    
    /**
     * 
     * @return readModel
     */
    private String getReadModel(){
        return readModel;
    }
    
    /**
     * 
     * @return mixedReads
     */
    private String getMixedReads(){
        return mixedReads;
    }

    /**
     * 
     * @param strMixedReads set
     */
    private void setMixedReads(String strMixedReads) {
        this.mixedReads = strMixedReads;
    }   
}
