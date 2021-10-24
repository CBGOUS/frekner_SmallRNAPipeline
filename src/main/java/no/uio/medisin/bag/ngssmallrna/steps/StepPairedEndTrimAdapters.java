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
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Performs adapter trimming of paired end reads
 * 
 * Not sure how well this is tested
 * 
 * @author joey
 */
public class StepPairedEndTrimAdapters extends NGSStep{

    static Logger logger = LogManager.getLogger();

    public static final NGSStepSubclass STEP_SUBCLASS           = NGSStepSubclass.DATABLE;
    
    public static final String STEP_ID_STRING = "PairedReadsAdaptersTrim";
    
    private static final String ID_SOFTWARE = "adapterSoftware";
    private static final String ID_ADAPTER_FILES = "adapterFile";
    
    private static final String ID_MISMATCHES ="noOfMismatches";
    private static final String ID_MIN_ALIGN_SCORE="noOfMismatches";
    private static final String ID_MIN_AVGRED_QUAL ="minAvgReadQual";
    
    private static final String ID_THREADS="noOfThreads";
    
    private static final String INFILES_EXTENSION = ".fastq";
    private static final String OUTFILES_EXTENSION = ".trim.paired.fastq";
    private static final String OUTFILES_EXTENSION_UNPAIRED = ".trim.unpaired.fastq";
    
    private static final String        RAW_INPUT_EXTENSION      = ".fastq.gz";
    
    private String pathToHiSAT="";
    private String adapterFile="";
    
    //input filename for forward sequence and reverse sequence files 
    private String fastqFile1 = "";
    private String fastqFile2 = "";
    
    
    
    //TODO this parameters are from single end trim class
    //     make sure they are using the same ones.
    private int noOfThreads=1;
    private int noOfMismatches=2;
    private int minALignScore=7;
    private int minAvgReadQuality=30;




    public StepPairedEndTrimAdapters(){
        classSubtype = NGSStepSubclass.DATABLE;
    }
    
    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepPairedEndTrimAdapters(InputDataForStep sid){
        classSubtype = NGSStepSubclass.DATABLE;
        stepInputData = sid;
    }

    
    
    @Override
    public String shortStepDescription(){
      return "Performs adapter trimming of paired end reads";
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "Performs adapter trimming of paired end reads.\n"
              + "The step requires the Trimmomatic software to be installed.\n"
              + "The location can be specified in the YAML configuration file\n"
              + "together with other parameters\n";
    }
    
    
    @Override
    public void parseStepParameters() throws Exception{
      
    }
    
    
    
    @Override
    public void verifyInputData() throws IOException {
        logger.info(STEP_ID_STRING + ": verify paired end trimming input data");
        
        this.setPaths();
        
        //check Trimming software is present
        if (new File(this.getPathToHiSAT()).exists() == false) {
            logger.error(STEP_ID_STRING + ": Adapter Trimming software is not found at location < " +this.getPathToHiSAT() + " >");
            throw new IOException(STEP_ID_STRING + ": Adapter Trimming software is not found at location < " +this.getPathToHiSAT() + " >");  
        }
        
        //check adapter files
        if(new File(this.getAdapterFile()).exists() == false){
            logger.error(STEP_ID_STRING + ": Adapter sequence files are not found at location < " +this.getAdapterFile() + " >");
            throw new IOException(STEP_ID_STRING + ": Adapter sequence files are not found at location < " +this.getAdapterFile() + " >"); 
        }
        
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        while(itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            fastqFile1 = inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, INFILES_EXTENSION );
            fastqFile2 = inFolder + FILESEPARATOR + sampleData.getFastqFile2().replace(RAW_INPUT_EXTENSION, INFILES_EXTENSION);
            
            //check fastq file 1 and 2 
            if(fastqFile1 == null){
                logger.error(STEP_ID_STRING + " : no Fastq file 1 forward specified");
                throw new IOException(STEP_ID_STRING + " : no fastq file 1 specified");
            }
            else if (fastqFile2 == null){
                logger.error(STEP_ID_STRING + ": no fastq file 2 specified");
                throw new IOException(STEP_ID_STRING + ": no fastq file 2 specified!");
            }
            
            //check files exist
            if((new File(this.cleanPath(fastqFile1))).exists() == false ){
                logger.error("AdapterTrimming: fastq file 1 <" + fastqFile1 +"> does not exist");
                throw new IOException("AdapterTrimming: fastq file 1 <" + fastqFile1 +"> does not exist");
            }
            else if((new File(this.cleanPath(fastqFile2))).exists() == false){
                logger.error("AdapterTrimming: fastq file 2 <" + fastqFile2 +"> does not exist");
                throw new IOException("AdapterTrimming: fastq file 2 <" + fastqFile2 +"> does not exist");
            }
            
            //check file extension
            if(fastqFile1.toUpperCase().endsWith(INFILES_EXTENSION.toUpperCase()) == false){
                logger.error(STEP_ID_STRING + ": incorrect file extension from input file 1 <" + fastqFile1 + ">." 
                +"should have <" + INFILES_EXTENSION + "> as extension");
                
                throw new IOException(STEP_ID_STRING + ": incorrect file extension from input file 1 <" + fastqFile1 + ">." 
                +"should have <" + INFILES_EXTENSION + "> as extension");
            } 
            else if(fastqFile2.toUpperCase().endsWith(INFILES_EXTENSION.toUpperCase()) == false){
                logger.error(STEP_ID_STRING + ": incorrect file extension from input file 2 <" + fastqFile1 + ">." 
                +"should have <" + INFILES_EXTENSION + "> as extension");
                
                throw new IOException(STEP_ID_STRING + ": incorrect file extension from input file 2 <" + fastqFile1 + ">." 
                +"should have <" + INFILES_EXTENSION + "> as extension");
            } 
        }
    }

    @Override
    public void verifyOutputData() throws  IOException{
        String outputFile1Paired = fastqFile1.replace(INFILES_EXTENSION, OUTFILES_EXTENSION);
        String outputFile1Unpaired = fastqFile1.replace(INFILES_EXTENSION, OUTFILES_EXTENSION_UNPAIRED);
        
        String outputFile2Paired = fastqFile2.replace(INFILES_EXTENSION, OUTFILES_EXTENSION);
        String outputFile2Unpaired = fastqFile2.replace(INFILES_EXTENSION, OUTFILES_EXTENSION_UNPAIRED);


        //check output files exist
        if ((new File(this.cleanPath(outputFile1Paired))).exists() == false) {
            logger.error("AdapterTrimming: output forward paired sequence file <" + outputFile1Paired + "> does not exist");
            throw new IOException("AdapterTrimming: output forward paired sequence file <" + outputFile1Paired + "> does not exist");
        }
        
        if ((new File(this.cleanPath(outputFile1Unpaired))).exists() == false) {
            logger.error("AdapterTrimming: output forward unpaired sequence file <" + outputFile1Unpaired + "> does not exist");
            throw new IOException("AdapterTrimming: output forward unpaired sequence file <" + outputFile1Unpaired + "> does not exist");
        } 
        
        if ((new File(this.cleanPath(outputFile2Paired))).exists() == false) {
            logger.error("AdapterTrimming: output reverse paired sequence file <" + outputFile2Paired + "> does not exist");
            throw new IOException("AdapterTrimming: output reverse paired sequence file <" + outputFile2Paired + "> does not exist");
        } 
        
        if ((new File(this.cleanPath(outputFile2Unpaired))).exists() == false) {
            logger.error("AdapterTrimming: output reverse unpaired sequence file <" + outputFile2Unpaired + "> does not exist");
            throw new IOException("AdapterTrimming: output reverse unpaired sequence file <" + outputFile2Unpaired + "> does not exist");
        }        
    }

    
    
    
    
    
    @Override
    public void parseConfigurationData(HashMap configData) throws Exception {
        logger.info(STEP_ID_STRING + ": verify configuration data ");

        //checking ID_SOFTWARE
        if(configData.get(ID_SOFTWARE) == null) {
            logger.error("<" + configData.get(ID_SOFTWARE) + "> : missing definition in configuration file");
            throw new NullPointerException("<" + configData.get(ID_SOFTWARE) + "> : missing definition in configuration file");
        }
        
        //checking ID_ADAPTER_FILES
        if(configData.get(ID_ADAPTER_FILES) == null){
            logger.error("<" + configData.get(ID_ADAPTER_FILES) + "> : missing definition in configuration file");
            throw new NullPointerException("<" + configData.get(ID_ADAPTER_FILES) + "> : missing definition in configuration file");
        }
        
        //checking ID_MISMATCHES
        if(configData.get(ID_MISMATCHES) == null){
            logger.error("<" + configData.get(ID_MISMATCHES) + "> : missing definition in configuration file");
            throw new NullPointerException("<" + configData.get(ID_MISMATCHES) + "> : missing definition in configuration file");
        }
        
        //checking ID_MIN_ALIGN_SCORE
        if(configData.get(ID_MIN_ALIGN_SCORE) == null){
            logger.error("<" + configData.get(ID_MIN_ALIGN_SCORE) + "> : missing definition in configuration file");
            throw new NullPointerException("<" + configData.get(ID_MIN_ALIGN_SCORE) + "> : missing definition in configuration file");
        }
        
        //checking ID_MIN_AVGRED_QUAL
        if(configData.get(ID_MIN_AVGRED_QUAL) == null){
            logger.error("<" + configData.get(ID_MIN_AVGRED_QUAL) + "> : missing definition in configuration file");
            throw new NullPointerException("<" + configData.get(ID_MIN_AVGRED_QUAL) + "> : missing definition in configuration file");
        }
        
        //checking ID_THREADS
        if(configData.get(ID_THREADS) == null){
            logger.error("<" + configData.get(ID_THREADS) + "> : missing definition in configuration file");
            throw new NullPointerException("<" + configData.get(ID_THREADS) + "> : missing definition in configuration file");
        }
        
        /**
         * this is different from single end class by simon.
         * checking ID_MISMATCHES
         * first is to check if it is a positive number,
         * if not, then won't call the method setNoOfMismatches.
         */
        if(((Integer)configData.get(ID_MISMATCHES)) <= 0){
            logger.error(ID_MISMATCHES + "<" + configData.get(ID_MISMATCHES) + "> must be positive");
            throw new IllegalArgumentException(ID_MISMATCHES + "<" + configData.get(ID_MISMATCHES) + "> must be positive");
        }
        else{
            try {
                this.setNoOfMismatches((Integer)configData.get(ID_MISMATCHES));
            } 
            catch (NumberFormatException exNm) {
                logger.error(ID_MISMATCHES + "<" + configData.get(ID_MISMATCHES) + "> is not an integer" );
                throw new NumberFormatException(ID_MISMATCHES + "<" + configData.get(ID_MISMATCHES) + "> is not an integer" );
            }
        }
        
        /**
         * checking ID_MIN_ALIGN_SCORE
         */
        try{
            this.setMinAlignScore((Integer) configData.get(ID_MIN_ALIGN_SCORE));
        }
        catch(NumberFormatException exNm){
            logger.error(ID_MIN_ALIGN_SCORE + " <" + configData.get(ID_MIN_ALIGN_SCORE) + "> is not an integer");
            throw new NumberFormatException(ID_MIN_ALIGN_SCORE + " <" + configData.get(ID_MIN_ALIGN_SCORE) + "> is not an integer");
        }        
        if (this.getMinAlignScore() <= 0){
            logger.error(ID_MIN_ALIGN_SCORE + " <" + configData.get(ID_MIN_ALIGN_SCORE) + "> must be positive");
            throw new IllegalArgumentException(ID_MIN_ALIGN_SCORE + " <" + configData.get(ID_MIN_ALIGN_SCORE) + "> must be positive");
        }    
        
        if(((Integer)configData.get(ID_MIN_ALIGN_SCORE)) <= 0){
            logger.error(ID_MIN_ALIGN_SCORE + "<" + configData.get(ID_MIN_ALIGN_SCORE) + "> must be positive");
            throw new IllegalArgumentException(ID_MIN_ALIGN_SCORE + "<" + configData.get(ID_MIN_ALIGN_SCORE) + "> must be positive");
        }
        else{
            try {
                this.setMinAlignScore((Integer)configData.get(ID_MIN_ALIGN_SCORE));
            } 
            catch (NumberFormatException exNm) {
                logger.error(ID_MIN_ALIGN_SCORE + "<" + configData.get(ID_MIN_ALIGN_SCORE) + "> is not an integer" );
                throw new NumberFormatException(ID_MIN_ALIGN_SCORE + "<" + configData.get(ID_MIN_ALIGN_SCORE) + "> is not an integer" );
            }
        }
        
        /**
         * checking ID_MIN_AVGRED_QUAL
         */
        if(((Integer)configData.get(ID_MIN_AVGRED_QUAL)) <= 0){
            logger.error(ID_MIN_AVGRED_QUAL + "<" + configData.get(ID_MIN_AVGRED_QUAL) + "> must be positive");
            throw new IllegalArgumentException(ID_MIN_AVGRED_QUAL + "<" + configData.get(ID_MIN_AVGRED_QUAL) + "> must be positive");
        }
        else{
            try {
                this.setMinAvgReadQuality((Integer)configData.get(ID_MIN_AVGRED_QUAL));
            } 
            catch (NumberFormatException exNm) {
                logger.error(ID_MIN_AVGRED_QUAL + "<" + configData.get(ID_MIN_AVGRED_QUAL) + "> is not an integer" );
                throw new NumberFormatException(ID_MIN_AVGRED_QUAL + "<" + configData.get(ID_MIN_AVGRED_QUAL) + "> is not an integer" );
            }
        }
        
        /**
         * checking ID_THREADS
         */
        if(((Integer)configData.get(ID_THREADS)) <= 0){
            logger.error(ID_THREADS + "<" + configData.get(ID_THREADS) + "> must be positive");
            throw new IllegalArgumentException(ID_THREADS + "<" + configData.get(ID_THREADS) + "> must be positive");
        }
        else{
            try {
                this.setNoOfThreads((Integer)configData.get(ID_THREADS));
            } 
            catch (NumberFormatException exNm) {
                logger.error(ID_THREADS + "<" + configData.get(ID_THREADS) + "> is not an integer" );
                throw new NumberFormatException(ID_THREADS + "<" + configData.get(ID_THREADS) + "> is not an integer" );
            }
        }
        
        
        this.setTrimSoftware((String)configData.get(ID_SOFTWARE));
        this.setAdapterFile((String) configData.get(ID_ADAPTER_FILES));
        logger.info("passed");
    }

    @Override
    public HashMap generateExampleConfigurationData() {
        logger.info(STEP_ID_STRING + ": generate example configuration data");
        HashMap configData = new HashMap();
        
        configData.put(StepPairedEndTrimAdapters.ID_SOFTWARE, "/home/joey/oslo/softwaresNGS/trimmomatic-0.36.jar");
        configData.put(StepPairedEndTrimAdapters.ID_ADAPTER_FILES, "/home/joey/oslo/softwareNGS/adapters/TruSeq2-PE.fa");
        configData.put(StepPairedEndTrimAdapters.ID_MISMATCHES, 2);
        configData.put(StepPairedEndTrimAdapters.ID_MIN_ALIGN_SCORE, 7);
        configData.put(StepPairedEndTrimAdapters.ID_MIN_AVGRED_QUAL,  30);
        configData.put(StepPairedEndTrimAdapters.ID_THREADS, 1);
        
        return configData;
    }

    @Override
    public void execute() throws IOException {
        logger.info(STEP_ID_STRING + ": execute step");
        
        String cmdTrimAdapters = "";
        
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        while (itSD.hasNext()){
            try{
                SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
                ArrayList<String> cmd = new ArrayList<>();
                
                cmd.add("java -jar");
                cmd.add(this.getPathToHiSAT());
                cmd.add("PE");
                cmd.add("-phred64 ");
                cmd.add("-threads " + this.getNoOfThreads());
                //todo
                cmd.add(inFolder + FILESEPARATOR + sampleData.getFastqFile1());
                cmd.add(inFolder + FILESEPARATOR + sampleData.getFastqFile2());
                
                Boolean f = new File(outFolder).mkdir();
                if(f)logger.info("created output folder <" + outFolder + "> for results");
                
                cmd.add(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, OUTFILES_EXTENSION));
                cmd.add(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, OUTFILES_EXTENSION_UNPAIRED));
                cmd.add(outFolder + FILESEPARATOR + sampleData.getFastqFile2().replace(RAW_INPUT_EXTENSION, OUTFILES_EXTENSION));
                cmd.add(outFolder + FILESEPARATOR + sampleData.getFastqFile2().replace(RAW_INPUT_EXTENSION, OUTFILES_EXTENSION_UNPAIRED));
                
                cmd.add("ILLUMINACLIP:" + this.getAdapterFile() + ": " + this.getNoOfMismatches() + "30" + ":" + this.getMinAlignScore());
                
                cmdTrimAdapters = this.cleanPath(StringUtils.join(cmd, " "));
                logger.info("Adapter Trim command:\t" + cmdTrimAdapters);
                
                Runtime runtime = Runtime.getRuntime();
                Process process = runtime.exec(cmdTrimAdapters);
                
                BufferedReader brStdIn = new BufferedReader(new InputStreamReader(process.getInputStream()));
                BufferedReader brStdErr = new BufferedReader(new InputStreamReader(process.getErrorStream()));
                
                String line = null;
                logger.info("<OUTPUT>");
                while ((line = brStdIn.readLine()) != null){
                    logger.info(line);
                }
                logger.info("</output>");
                
                logger.info("<ERROR>");
                while((line = brStdErr.readLine()) != null) {
                    logger.info(line);
                }
                logger.info("</ERROR>");
                
                int exitVal = process.waitFor();
                logger.info("Process exitValue: " + exitVal);
                
                brStdIn.close();
                brStdErr.close();
                
            } 
            catch (IOException|InterruptedException ex) {
                logger.error("error executing AdapterTrimming command\n" + ex.toString());
                throw new IOException(STEP_ID_STRING + "error executing AdapterTrimming command " + cmdTrimAdapters);
            }
        }
        logger.info(STEP_ID_STRING + ": completed");
    }

    
    
    
    @Override
    public NGSStepSubclass getStepSubclass(){
        return STEP_SUBCLASS;
    }
    
    
    
    
    /**
     * @return pathToHiSAT
     */
    private String getPathToHiSAT() {
        return pathToHiSAT; 
    }
    
    /**
     * 
     * @param trimSoftware 
     */
    private void setTrimSoftware(String trimSoftware){
        this.pathToHiSAT = trimSoftware;
    }

  
    /**
     * 
     * @return adapterFile
     */
    private String getAdapterFile() {
        return adapterFile;
    }
    
    /**
     * 
     * @param adapterFile set 
     */
    private void setAdapterFile(String adapterFiles){
        this.adapterFile = adapterFiles;
    }

    /**
     * 
     * @param noOfMismatches set
     */
    private void setNoOfMismatches(int noOfMismatches) {
        this.noOfMismatches = noOfMismatches;
    }

    /**
     * 
     * @param minAlignScore  set
     */
    private void setMinAlignScore(int minAlignScore) {
        this.minALignScore = minAlignScore;
    }

    /**
     * 
     * @param minAvgReadQuality set
     */
    private void setMinAvgReadQuality(int minAvgReadQuality) {
        this.minAvgReadQuality = minAvgReadQuality;
    }

    /**
     * 
     * @param noOfThreads set
     */
    private void setNoOfThreads(int noOfThreads) {
        this.noOfThreads = noOfThreads; 
    }
    
    /**
     * 
     * @return noOfThreads
     */
    private int getNoOfThreads() {
        return noOfThreads;
    }

    /**
     * 
     * @return noOfMismatches
     */
    private int getNoOfMismatches() {
        return noOfMismatches;
    }

    /**
     * 
     * @return minALignScore
     */
    private int getMinAlignScore() {
        return minALignScore;
    }


   
}
