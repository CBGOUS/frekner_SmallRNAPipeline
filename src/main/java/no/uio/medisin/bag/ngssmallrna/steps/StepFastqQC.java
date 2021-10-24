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
//import javax.xml.bind.ValidationException;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import static no.uio.medisin.bag.ngssmallrna.steps.NGSStep.FILESEPARATOR;
import org.apache.commons.cli.ParseException;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;

import org.apache.logging.log4j.Logger;




/**
 *  FastQC Step
 *  perform QC of FASTQ file list.
 * 
 *   Input is a FASTQ file
 *   Output is a quality report
 * 
 * 
 * @author sr/joey
 */

public class StepFastqQC extends NGSStep {
    
    static Logger                       logger = LogManager.getLogger();
    
    public static final NGSStepSubclass STEP_SUBCLASS           = NGSStepSubclass.DATABLE;

    public static final String          STEP_ID_STRING          = "FastqQC";
    public static final String          ID_SOFTWARE             = "fastQCSoftware";    
    
    private static final String         INFILE_EXTENSION        = ".fastq";
    private static final String         RAW_INPUT_EXTENSION     = ".fastq.gz";
    
    private String                      fastqcSoftware          = "";

    

    public StepFastqQC(){
        classSubtype = NGSStepSubclass.DATABLE;
    }
    
    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepFastqQC(InputDataForStep sid){
        classSubtype = NGSStepSubclass.DATABLE;
        stepInputData = sid;
    }
    
    

    @Override
    public String shortStepDescription(){
      return "Performs QC analysis of each specified FASTQ file";
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "unzip input FASTQ files from gzip format.\n"
              + "The step requires the FASTQC software to be installed.\n"
              + "The location can be specified in the YAML configuration file\n"
              + "The step uses FASTQC in single mode, rather than batch analysis "
              + "so that a separate report is generated for each each file\n";
    }
    
    
    @Override
    public void parseStepParameters() throws Exception{
      
    }
    
    

    
    /**
     * This parses out the run parameters for this step
     * 
     * @param configData
     * @throws ValidationException 
     */
    @Override
    public void parseConfigurationData(HashMap configData) throws ParseException{

        logger.info(STEP_ID_STRING + ": verify configuration data");
        
        if(configData.get(ID_SOFTWARE)==null) {
            logger.error("<" + configData.get(ID_SOFTWARE) + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + configData.get(ID_SOFTWARE) + "> : Missing Definition in Configuration File");
        }

        this.setFastQCSoftware((String) configData.get(ID_SOFTWARE));
 
        logger.info("passed");
    }
    
    
    
    /**
     * fastqc specified file types
     * fastqc can take both .fastq and .fastq.gz format 
     * but here we only use .fastq.gz, before unzip reads files.
     * @throws IOException 
     */
    @Override
    public void execute() throws IOException{
        /*
            sample fastqcSoftware command       
            fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam]  
         */
        logger.info(STEP_ID_STRING + ": execute step");

        Boolean fA = new File(outFolder).mkdir();
        if (fA) {
            logger.info("created output folder <" + outFolder + "> for results");
        }

        String fastqFile1in = "";
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        while (itSD.hasNext()) {
            String fastqcCmdString = "";
            try {
                SampleDataEntry sampleData = (SampleDataEntry) itSD.next();

                fastqFile1in = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION));

                ArrayList<String> cmd = new ArrayList<>();
                cmd.add(this.getFastQCSoftware());
                cmd.add("-o");
                cmd.add(outFolder);
                cmd.add(fastqFile1in);
                if (sampleData.getFastqFile2() != null) {
                    cmd.add(inFolder + FILESEPARATOR + sampleData.getFastqFile2().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION));
                }
                fastqcCmdString = this.cleanPath(StringUtils.join(cmd, " "));
                logger.info("fastqc command:\t" + fastqcCmdString);
                Runtime rt = Runtime.getRuntime();
                Process proc = rt.exec(fastqcCmdString);
                BufferedReader brAStdin = new BufferedReader(new InputStreamReader(proc.getInputStream()));
                BufferedReader brAStdErr = new BufferedReader(new InputStreamReader(proc.getErrorStream()));

                String line = null;
                logger.info("<OUTPUT>");
                while ((line = brAStdin.readLine()) != null) {
                    logger.info(line);
                }
                logger.info("</OUTPUT>");

                logger.info("<ERROR>");
                while ((line = brAStdErr.readLine()) != null) {
                    logger.info(line);
                }
                logger.info("</ERROR>");

                int exitVal = proc.waitFor();
                logger.info("Process exitValue: " + exitVal);

                brAStdin.close();
                brAStdErr.close();
                logger.info("FastQC done");

            } catch (Exception ex) {
                logger.error("error executing fastqc command");
                throw new IOException("error executing fastqc command");
            }
        }
        
        logger.info(STEP_ID_STRING + ": completed");
    }
    
    
    
            
    /**
     * this should be called prior to executing the step.
     * check fastqc software exists and input files are available
     * 
     * @throws IOException
     */
    @Override
    public void verifyInputData() throws IOException{
        
        logger.info("verify input data");        
        this.setPaths();
        
        if(new File(this.getFastQCSoftware()).exists() == false){
            logger.error("fastqc software not found at location < " + this.getFastQCSoftware() +">");
            throw new IOException("fastqc software not found at location < " + this.getFastQCSoftware() +">");
        }
                
        // check the data files
        String fastqFile1in = "";
        String fastqFile1out = "";
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            
            //Fastq 1
            if (sampleData.getFastqFile1()==null) {
                logger.error("no Fastq1 file specified");
                throw new IOException("no Fastq1 file specified");
            }
            
            fastqFile1in = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION));
            if (new File(fastqFile1in).exists()==false && new File(fastqFile1in.replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION)).exists()==false){
                logger.error(STEP_ID_STRING + ": fastq 1 file 1 <" + fastqFile1in + "> & <" + fastqFile1out + "> do not exist");
                throw new IOException(STEP_ID_STRING + ": fastq file 1 <" + fastqFile1in + "> & <" + fastqFile1out + "> does not exist");
            }
            if (fastqFile1in.replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION).toUpperCase().endsWith(INFILE_EXTENSION.toUpperCase())==false)
            {
                logger.info(STEP_ID_STRING + ": incorrect file extension for input file <" 
                  + fastqFile1in.replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION) + ">. should have <" + INFILE_EXTENSION + "> as extension");
                throw new IOException(STEP_ID_STRING + ": incorrect file extension for input file <" 
                  + fastqFile1in.replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION) + ">.  \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
            }
 
            //Fastq 2
            if (sampleData.getFastqFile2()==null) continue;
            String fastqFile2in = "";
            String fastqFile2out = "";
            fastqFile2in = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile2().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION));

            
            if ((new File(fastqFile2in)).exists()==false && (new File(fastqFile2in)).exists()==false){
                logger.error(STEP_ID_STRING + ": fastq 2 Files <" + fastqFile2in + "> & <" + fastqFile2out + "> do not exist");
                throw new IOException(STEP_ID_STRING + ": fastq 2 Files <" + fastqFile2in + "> & <" + fastqFile2out + "> do not exist");
            }
            if (fastqFile2in.toUpperCase().endsWith(INFILE_EXTENSION.toUpperCase())==false)
            {
                logger.error(STEP_ID_STRING + ": incorrect file extension for fastq file 2 <" 
                  + fastqFile2in + ">. should have <" + INFILE_EXTENSION + "> as extension");
                throw new IOException(STEP_ID_STRING + ": incorrect file extension for fastq file 2 <" 
                  + fastqFile2in + ">. \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
            }                       
        }
    }
    
    
    
    @Override
    public void verifyOutputData(){
        logger.info("no output verification required");
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
        
        HashMap<String, Object> configData = new HashMap();
        
        configData.put(ID_SOFTWARE, "/usr/local/fastqc");
        
        return configData;
        
    }

    
    
    
    
    
    /**
     * @return the fastqcSoftware
     */
    public String getFastQCSoftware() {
        return fastqcSoftware;
    }

    /**
     * @param fastqcSoftware the fastqcSoftware to set
     */
    public void setFastQCSoftware(String fastqcSoftware) {
        this.fastqcSoftware = fastqcSoftware;
    }

    
}
