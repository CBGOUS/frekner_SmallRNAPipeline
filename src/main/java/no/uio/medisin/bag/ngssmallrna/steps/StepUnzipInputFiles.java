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
 *  Adapter Trimming Step
 *  Unzip FASTQ files using pigz.
 * 
 *    Input is a zipped FASTQ file
 *    Output is a unzipped FASTQ file
 * 
 *  06/10/2016
 *  updated to use the checkParameters() method to check parameter value is valid
 * 
 * @author sr
 */

public class StepUnzipInputFiles extends NGSStep{
    
    static Logger                   logger                          = LogManager.getLogger();
    
    public static final NGSStepSubclass STEP_SUBCLASS           = NGSStepSubclass.DATABLE;
    
    public static final String          STEP_ID_STRING          = "UnzipInputFiles";
    private static final String         ID_SOFTWARE             = "unzipSoftware";    
    private static final String         ID_THREADS              = "noOfThreads";
    
    private static final String         INFILE_EXTENSION        = ".fastq.gz";
    private static final String         OUTFILE_EXTENSION       = ".fastq";
    private static final String         RAW_INPUT_EXTENSION     = ".fastq.gz";
    
    private int                         noOfThreads             = 4;
    private String                      unzipSoftware           = "";

    

    public StepUnzipInputFiles(){
       classSubtype = NGSStepSubclass.DATABLE;
    }
    
    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepUnzipInputFiles(InputDataForStep sid){
       stepInputData = sid;
       classSubtype = NGSStepSubclass.DATABLE;
    }
    
    
    @Override
    public String shortStepDescription(){
      return "unzip input FASTQ files from gzip format";
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "unzip input FASTQ files from gzip format.\n"
              + "The step requires the PIGZ software to be installed.\n"
              + "The location can be specified in the YAML configuration file\n\n";
    }
    
    @Override
    public void parseStepParameters() throws Exception{
      
    }
    
    
    
    /**
     * Parse out the HashMap containing the run parameters for this step
     * 
     * @param configData
     * @throws Exception 
     */
    @Override
    public void parseConfigurationData(HashMap configData) throws Exception{

        logger.info(STEP_ID_STRING + ": verify configuration data");
        
        if(configData.get(ID_SOFTWARE)==null) {
            logger.error("<" + ID_SOFTWARE + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_SOFTWARE + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_THREADS)==null) {
            logger.error("<" + ID_THREADS + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_THREADS + "> : Missing Definition in Configuration File");
        }
        
        String chk;

        chk = checkParameter("Integer", ID_THREADS, Integer.toString((Integer)configData.get(ID_THREADS)), "0", "12", logger);
        if(chk!=null)
            this.setNoOfThreads((Integer)configData.get(ID_THREADS));
            

        this.setUnzipSoftware((String) configData.get(ID_SOFTWARE));
 
        logger.info("passed");
    }
   
    /**
     * unzip specified file types
     * 
     * @throws IOException 
     */
    @Override
    public void execute() throws IOException{
        /*
            sample unzip command       
            pigz -p 4 -d /data/ngsdata/project1/sra_data.fastq.gz     
        */
        logger.info(STEP_ID_STRING + ": execute step");        
        
        String fastqFile1in = "";
        String fastqFile1out = "";
        String fastqFile2in = "";
        String fastqFile2out = "";
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        while (itSD.hasNext()){
            try{
                SampleDataEntry sampleData = (SampleDataEntry)itSD.next();                
                
                fastqFile1out = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(INFILE_EXTENSION, OUTFILE_EXTENSION));
                
                if (new File(sampleData.getFastqFile1()).exists()){
                    fastqFile1in = this.cleanPath(sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION));
                }
                else{
                    fastqFile1in = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION));
               
                } 
                if(new File(fastqFile1out).exists()==false){
                                
                    ArrayList<String> cmd1 = new ArrayList<>();
                    cmd1.add(this.getUnzipSoftware());
                    cmd1.add("-d");
                    cmd1.add(fastqFile1in);
                    cmd1.add("-p " + this.getNoOfThreads());

                    String cmdUnzip = StringUtils.join(cmd1, " ");
                    cmdUnzip = cmdUnzip.replace(FILESEPARATOR + FILESEPARATOR, FILESEPARATOR);
                    logger.info(STEP_ID_STRING + " fastq1 unzip command:\t" + cmdUnzip);

                    Runtime rt1 = Runtime.getRuntime();
                    Process proc1 = rt1.exec(cmdUnzip);
                    BufferedReader brStdin1  = new BufferedReader(new InputStreamReader(proc1.getInputStream()));
                    BufferedReader brStdErr1 = new BufferedReader(new InputStreamReader(proc1.getErrorStream()));

                        String line1 = null;
                        logger.info("<OUTPUT>");
                        while ( (line1 = brStdin1.readLine()) != null)
                            logger.info(line1);
                        logger.info("</OUTPUT>");

                        logger.info("<ERROR>");
                        while ( (line1 = brStdErr1.readLine()) != null)
                            logger.info(line1);
                        logger.info("</ERROR>");                

                    int exitVal = proc1.waitFor();            
                    logger.info("Process1 exitValue: " + exitVal);   

                    brStdin1.close();
                    brStdErr1.close();
                }
                else{               
                    logger.info("fastq file 1 <" + this.cleanPath(fastqFile1out) + "> exists. Skipping");
                }
                
                if (sampleData.getFastqFile2() == null) continue;
                fastqFile2in = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile2().replace(".fastq", INFILE_EXTENSION));
                fastqFile2out = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile2());
                if(new File(fastqFile2out).exists()==false){
            
                    ArrayList<String> cmd2 = new ArrayList<>();
                    cmd2.add(this.getUnzipSoftware());
                    cmd2.add("-d");
                    cmd2.add(fastqFile2in);
                    cmd2.add("-p " + this.getNoOfThreads());

                    String cmdUnzip2 = StringUtils.join(cmd2, " ");
                    cmdUnzip2 = cmdUnzip2.replace(FILESEPARATOR + FILESEPARATOR, FILESEPARATOR);
                    logger.info(STEP_ID_STRING + " fastq2 unzip command:\t" + cmdUnzip2);

                    Runtime rt2 = Runtime.getRuntime();
                    Process proc2 = rt2.exec(cmdUnzip2);
                    BufferedReader brStdin2  = new BufferedReader(new InputStreamReader(proc2.getInputStream()));
                    BufferedReader brStdErr2 = new BufferedReader(new InputStreamReader(proc2.getErrorStream()));

                        String line2 = null;
                        logger.info("<OUTPUT>");
                        while ( (line2 = brStdin2.readLine()) != null)
                            logger.info(line2);
                        logger.info("</OUTPUT>");

                        logger.info("<ERROR>");
                        while ( (line2 = brStdErr2.readLine()) != null)
                            logger.info(line2);
                        logger.info("</ERROR>");                

                    int exitVal2 = proc2.waitFor();            
                    logger.info("Process2 exitValue: " + exitVal2);   

                    brStdin2.close();
                    brStdErr2.close();
                }
                else{               
                    logger.info("fastq file 2 <" + this.cleanPath(fastqFile2out) + "> exists. Skipping");
                }
                                
            }
            catch(InterruptedException ex){
                logger.error("error executing pigz unzip command");
                throw new IOException("error executing pigz unzip command");
            }
        }
        
        logger.info(STEP_ID_STRING + ": completed");
    }
    
    
    
            
    /**
     * this should be called prior to executing the step.
     * check unzip software exists and input files are available
     * 
     * @throws IOException
     */
    @Override
    public void verifyInputData() throws IOException{
        
        logger.info("verify input data");        
        this.setPaths();
        
        if(new File(this.getUnzipSoftware()).exists() == false){
            logger.error("unzip software not found at location < " + this.getUnzipSoftware() +">");
            throw new IOException("unzip software not found at location < " + this.getUnzipSoftware() +">");
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
            
            fastqFile1out = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, OUTFILE_EXTENSION));
            fastqFile1in = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION));
            String fastqFile1 = inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION);
            if (new File(fastqFile1in).exists()==false && new File(fastqFile1out).exists()==false){
                logger.error(STEP_ID_STRING + ": fastq 1 Files <" + fastqFile1in + "> & <" + fastqFile1out + "> do not exist");
                throw new IOException(STEP_ID_STRING + ": fastq 1 Files <" + fastqFile1in + "> & <" + fastqFile1out + "> do not exist");
            }
            if (fastqFile1.toUpperCase().endsWith(INFILE_EXTENSION.toUpperCase())==false)
            {
                logger.info(STEP_ID_STRING + ": incorrect file extension for input file <" 
                  + fastqFile1 + ">. should have <" + INFILE_EXTENSION + "> as extension");
                throw new IOException(STEP_ID_STRING + ": incorrect file extension for input file <" 
                  + fastqFile1 + ">.  \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
            }
            
            
            //Fastq 2
            if (sampleData.getFastqFile2()==null) continue;
            String fastqFile2in = "";
            String fastqFile2out = "";
            fastqFile2in = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile2().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION));
            fastqFile2out = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile2());

            
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
        
        configData.put(ID_SOFTWARE, "/usr/local/pigz");
        configData.put(ID_THREADS, 4);
        
        return configData;
        
    }


    public String describeStep(){
      return "Unzip the input FASTQ files in preparation for subsequent analysis steps";
    }
    
    /**
     * @return the noOfThreads
     */
    public int getNoOfThreads() {
        return noOfThreads;
    }

    /**
     * @param noOfThreads the noOfThreads to set
     */
    public void setNoOfThreads(int noOfThreads) {
        this.noOfThreads = noOfThreads;
    }

    /**
     * @return the unzipSoftware
     */
    public String getUnzipSoftware() {
        return unzipSoftware;
    }

    /**
     * @param unzipSoftware the unzipSoftware to set
     */
    public void setUnzipSoftware(String unzipSoftware) {
        this.unzipSoftware = unzipSoftware;
    }

}
