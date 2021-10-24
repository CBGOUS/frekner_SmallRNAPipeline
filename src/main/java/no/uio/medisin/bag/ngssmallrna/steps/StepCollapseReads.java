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
 *  Collapse Read Step
 *  1. Convert FASTQ files in input list to an equivalent set of FASTA files.
 *  2. Count up duplicate reads and store this information in the header line 
 * 
 * @author sr
 */

public class StepCollapseReads extends NGSStep{
    
    static  Logger                      logger = LogManager.getLogger();
    
    public static final NGSStepSubclass STEP_SUBCLASS           = NGSStepSubclass.DATABLE;

    public static final String          STEP_ID_STRING          = "CollapseReads";
    private static final String         ID_Q2A_SOFTWARE         = "fastq2fastaLocation";   
    private static final String         ID_COLLAPSE_SOFTWARE    = "fastxCollapserLocation";   
    
    
    private static final String         INFILE_EXTENSION        = ".trim.fastq";
    private static final String         FASTA_OUTFILE_EXTENSION = ".trim.fasta";
    private static final String         CLP_OUTFILE_EXTENSION   = ".trim.clp.fasta";
    private static final String         RAW_INPUT_EXTENSION     = ".fastq.gz";
    
    private String                      fastq2fasta_software    = "";
    private String                      collapseFastaSoftware   = "";
    
    
    

    public StepCollapseReads(){
        classSubtype = NGSStepSubclass.DATABLE;
    }
    
    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepCollapseReads(InputDataForStep sid){
        classSubtype = NGSStepSubclass.DATABLE;
        stepInputData = sid;
    }
    
    
    @Override
    public String shortStepDescription(){
      return "Collapse FASTQ file to FASTA file";
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "Collapse FASTQ file to FASTA file.\n"
              + "The step requires the FASTX toolkit from CSHL to be installed.\n"
              + "The location can be specified in the YAML configuration file\n"
              + "Each FASTQ files is converted to an equivalent set of FASTA file.\n" 
              + "Duplicate reads are counted and the information is stored in the header line\n";
    }
    
    
    @Override
    public void parseStepParameters() throws Exception{
      
    }
    
    
    
    /**
     * This parses out the hashmap containing the run parameters for this step
     * 
     * @param configData
     * @throws Exception 
     */
    @Override
    public void parseConfigurationData(HashMap configData) throws Exception{

        logger.info(STEP_ID_STRING + ": verify configuration data");
    
        if(configData.get(ID_Q2A_SOFTWARE)==null) {
            logger.error("<" + configData.get(ID_Q2A_SOFTWARE) + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + configData.get(ID_Q2A_SOFTWARE) + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_COLLAPSE_SOFTWARE)==null) {
            logger.error("<" + configData.get(ID_COLLAPSE_SOFTWARE) + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + configData.get(ID_COLLAPSE_SOFTWARE) + "> : Missing Definition in Configuration File");
        }
        
        this.setFastq2fasta_software((String) configData.get(ID_Q2A_SOFTWARE));
        this.setCollapseFastaSoftware((String) configData.get(ID_COLLAPSE_SOFTWARE));
        
        logger.info("passed");
    }
    
    
    
    
    @Override
    public void execute() throws IOException{
        this.setPaths();
        
        String cmdFQ2FA = "";
        String cmdClp = "";
        
        String fastqInputFile = "";
        String fastaOutputFile = "";
        String clpOutputFile = ""; 
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            try{
                Boolean f = new File(outFolder).mkdir(); 
                if (f) logger.info("created output folder <" + outFolder + "> for results" );
 
                /*
                    fastq_to_fasta -i 1000.fastq -o 1000.fasta -Q33
                */                      
                ArrayList<String> cmdQ2A = new ArrayList<>();
                cmdQ2A.add(this.getFastq2fasta_software());

                fastqInputFile = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION));
                cmdQ2A.add("-i");
                cmdQ2A.add(fastqInputFile);

                cmdQ2A.add("-o");

                fastaOutputFile = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, FASTA_OUTFILE_EXTENSION);
                cmdQ2A.add(fastaOutputFile);
                
                cmdQ2A.add("-Q33");


                cmdFQ2FA = this.cleanPath(StringUtils.join(cmdQ2A, " "));
                logger.info("Fastq2Fasta command:\t" + cmdFQ2FA);
                Runtime rt = Runtime.getRuntime();
                Process proc = rt.exec(cmdFQ2FA);
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
                    System.out.println("Process exitValue: " + exitVal);            

                brStdin.close();
                brStdErr.close();
            }
            catch(IOException|InterruptedException exIE){
                logger.error("error executing Fastq2Fasta command\n");
                logger.error("CMD is " + cmdFQ2FA);
                throw new IOException(STEP_ID_STRING + ": error executing Fastq2Fasta command " + cmdFQ2FA);
            }
            
            
            
                /*
                    fastx_collapse -i 1000.fastq -o 1000.fasta -Q33
                */                      
            try{
                ArrayList<String> cmd2 = new ArrayList<>();
                cmd2.add(this.getCollapseFastaSoftware());

                String fastaInputFile = fastaOutputFile;
                cmd2.add("-i");
                cmd2.add(fastaInputFile);

                cmd2.add("-o");

                clpOutputFile = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, CLP_OUTFILE_EXTENSION);
                cmd2.add(clpOutputFile);
                

                cmdClp = this.cleanPath(StringUtils.join(cmd2, " "));
                logger.info("Collapse fasta command:\t" + cmdClp);

                Runtime rtClp = Runtime.getRuntime();
                Process procClp = rtClp.exec(cmdClp);
                BufferedReader brStdinClp  = new BufferedReader(new InputStreamReader(procClp.getInputStream()));
                BufferedReader brStdErrClp = new BufferedReader(new InputStreamReader(procClp.getErrorStream()));
                
                String line = null;
                logger.info("<OUTPUT>");
                while ( (line = brStdinClp.readLine()) != null)
                    logger.info(line);
                logger.info("</OUTPUT>");
                
                logger.info("<ERROR>");
                while ( (line = brStdErrClp.readLine()) != null)
                    logger.info(line);
                logger.info("</ERROR>");
                
                
                int exitValclp = procClp.waitFor();            
                System.out.println("Process exitValue: " + exitValclp);            
            
                brStdinClp.close();
                brStdErrClp.close();
            
            }
            catch(IOException|InterruptedException exIE){
                logger.error("error executing Collapse Fasta command\n");
                logger.error("CMD is " + cmdClp);
                throw new IOException(STEP_ID_STRING + ": error executing Collapse Fasta command " + cmdClp);
             }
        }
        
        logger.info(STEP_ID_STRING + ": completed");
        
    }
    

    
    /**
     * check data makes sense 
     */
    @Override
    public void verifyInputData() throws IOException, NullPointerException{
        
        logger.info(STEP_ID_STRING + ": verify input data");        
        this.setPaths();
                
        // does software exist?
        if(new File(this.getFastq2fasta_software()).exists() == false){
            logger.error(STEP_ID_STRING + ": fastq2fasta software not found at location < " + this.getFastq2fasta_software() +">");
            throw new IOException(STEP_ID_STRING + ": fastq2fasta software not found at location < " + this.getFastq2fasta_software() +">");
        }
        if(new File(this.getCollapseFastaSoftware()).exists() == false){
            logger.error(STEP_ID_STRING + ": fastq2fasta software not found at location < " + this.getCollapseFastaSoftware() +">");
            throw new IOException(STEP_ID_STRING + ": fastq2fasta software not found at location < " + this.getCollapseFastaSoftware() +">");
        }
        
        // check the data files
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            
            //Fastq 1
            if (sampleData.getFastqFile1()==null) {
                logger.error(STEP_ID_STRING + ": no Fastq1 file specified");
                throw new IOException(STEP_ID_STRING + ": no Fastq1 file specified");
            }
            String fastqFile1 = inFolder + NGSStep.FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION);
            
            if ((new File(fastqFile1)).exists()==false){
                logger.error(STEP_ID_STRING + ": fastq File1 <" 
                  + fastqFile1 + "> does not exist");
                throw new IOException(STEP_ID_STRING + " : fastq File1 <" 
                  + fastqFile1 + "> does not exist");
            }
            if (fastqFile1.toUpperCase().endsWith(INFILE_EXTENSION.toUpperCase())==false)
            {
                logger.error(STEP_ID_STRING + ": incorrect file extension for input file <" 
                  + fastqFile1 + ">.  \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
                throw new IOException(STEP_ID_STRING + ": incorrect file extension for input file <" 
                  + fastqFile1 + ">.  \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
            }
            
            
            //Fastq 2
            if (sampleData.getFastqFile2()==null) continue;
            String fastqFile2 = inFolder + NGSStep.FILESEPARATOR + sampleData.getFastqFile2().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION);
            
            if ((new File(fastqFile2)).exists()==false){
                logger.error(STEP_ID_STRING + " : fastq File2 <" 
                  + fastqFile2 + "> does not exist");
                throw new IOException(STEP_ID_STRING + " : fastq File2 <" 
                  + fastqFile2 + "> does not exist");
            }
            if (fastqFile2.toUpperCase().endsWith(INFILE_EXTENSION.toUpperCase())==false)
            {
                logger.error(STEP_ID_STRING + " : incorrect file extension for fastq file 2 <" 
                  + fastqFile2 + ">. \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
                throw new IOException(STEP_ID_STRING + " : incorrect file extension for fastq file 2 <" 
                  + fastqFile2 + ">. \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
            }
                        
            
        }

    }
    
    @Override
    public void verifyOutputData(){
        
    }

    
    @Override
    public NGSStepSubclass getStepSubclass(){
        return STEP_SUBCLASS;
    }
    
    

    
    /**
     * generate sample configuration data so the user can see what can be 
     * specified 
     * 
     * @return 
     */
    @Override
    public HashMap generateExampleConfigurationData() {

        logger.info(STEP_ID_STRING + ": generate example configuration data");
        
        HashMap configData = new HashMap();
        
        configData.put(ID_Q2A_SOFTWARE, "/usr/local/bin/fastq_to_fasta");
        configData.put(ID_COLLAPSE_SOFTWARE, "/usr/local/bin/fastx_collapser");
        
        return configData;
    }   
    
    
    
    
    
    /**
     * @return the fastq2fasta_software
     */
    public String getFastq2fasta_software() {
        return fastq2fasta_software;
    }

    /**
     * @param fastq2fasta_software the fastq2fasta_software to set
     */
    public void setFastq2fasta_software(String fastq2fasta_software) {
        this.fastq2fasta_software = fastq2fasta_software;
    }

    /**
     * @return the collapseFastaSoftware
     */
    public String getCollapseFastaSoftware() {
        return collapseFastaSoftware;
    }

    /**
     * @param collapseFastaSoftware the collapseFastaSoftware to set
     */
    public void setCollapseFastaSoftware(String collapseFastaSoftware) {
        this.collapseFastaSoftware = collapseFastaSoftware;
    }

    
}
