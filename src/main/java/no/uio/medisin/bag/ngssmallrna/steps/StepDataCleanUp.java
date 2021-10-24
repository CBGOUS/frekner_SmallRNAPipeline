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
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;

import org.apache.logging.log4j.Logger;


/**
 *  perform clean up
 *  this primarily involves zipping up files to save space
 *  because we don't know which steps were run or where the files are located
 *  we have to find files by searching the project folder
 *  thus, this will only work for linux command line, there is no attempt
 *  to make it work for Windows
 * 
 *  
 * 
 *   Input is a zipped FASTQ file
 *   Output is a unzipped FASTQ file
 * 
 * need to add the ability to define which types of files to compress. e.g., fastq, fasta, sam
 * 
 * @author sr
 */

public class StepDataCleanUp extends NGSStep{
    
    static Logger                       logger = LogManager.getLogger();

    public static final NGSStepSubclass STEP_SUBCLASS           = NGSStepSubclass.DATABLE;

    public  static final String         STEP_ID_STRING          = "DataCleanUp";
    private static final String         ID_SOFTWARE             = "unzipSoftware";    
    private static final String         ID_THREADS              = "noOfThreads";
    private static final String         ID_FILE_TYPES           = "fileTypes";
    
    
    
    private int                         noOfThreads             = 4;
    private String                      unzipSoftware           = "";
    private ArrayList<String>           fileTypes;
    
    
    

    public StepDataCleanUp(){
        classSubtype = NGSStepSubclass.DATABLE;
    }

    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepDataCleanUp(InputDataForStep sid){
        classSubtype = NGSStepSubclass.DATABLE;
        stepInputData = sid;
    }
    
    @Override
    public void parseStepParameters() throws Exception{
      
    }
    

    @Override
    public String shortStepDescription(){
      return "perform clean up";
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "perform clean up.\n\n"
              + "this primarily involves zipping up files to save space.\n "
              + "Because we don't know which steps were run or where the files are "
              + "located, we have to find files by searching the project folder. "
              + "Thus, this will only work for linux command line, there is no "
              + "attempt to make it work for Windows\n";
    }



    /**
     * This parses out the hashMap containing the run parameters for this step
     * 
     * @param configData
     * @throws Exception 
     */
    @Override
    public void parseConfigurationData(HashMap configData) throws Exception{

        logger.info(STEP_ID_STRING + ": verify configuration data");
        if(configData.get(ID_SOFTWARE)==null) {
            logger.error("<" + configData.get(ID_SOFTWARE) + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + configData.get(ID_SOFTWARE) + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_THREADS)==null) {
            logger.error("<" + configData.get(ID_THREADS) + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + configData.get(ID_THREADS) + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_FILE_TYPES)==null) {
            logger.error("<" + configData.get(ID_FILE_TYPES) + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + configData.get(ID_FILE_TYPES) + "> : Missing Definition in Configuration File");
        }
        
        
        try{
            this.setNoOfThreads((Integer)configData.get(ID_THREADS));
        }
        catch(NumberFormatException exNm){
            logger.error(ID_THREADS + " <" + configData.get(ID_THREADS) + "> is not an integer");
            throw new NumberFormatException(ID_THREADS + " <" + configData.get(ID_THREADS) + "> is not an integer");
        }
        
        if (this.getNoOfThreads() <= 0){
            logger.error(ID_THREADS + " <" + configData.get(ID_THREADS) + "> must be positive");
            throw new IllegalArgumentException(ID_THREADS + " <" + configData.get(ID_THREADS) + "> must be positive");
        }
        
        this.setUnzipSoftware((String) configData.get(ID_SOFTWARE));
        
        try{
            this.setFileTypes((ArrayList<String> )configData.get(ID_FILE_TYPES));
        }
        catch(Exception ex){
            logger.error("couldn't cast " + configData.get(ID_FILE_TYPES) + "value to ArrayList");
            throw new IOException("couldn't cast " + configData.get(ID_FILE_TYPES) + "value to ArrayList");
        }

        logger.info("passed");
    }
    
    
    
    
    @Override
    public void execute() throws IOException{
        
        logger.info(STEP_ID_STRING + ": execute step");        
        
        String projectRoot = getStepInputData().getProjectRoot() + FILESEPARATOR + getStepInputData().getProjectID();
        File directory = new File(projectRoot);
        String cmdZip = "";
        for(String fileType: fileTypes){
            String fileTypesArray[] = new String[fileTypes.size()];
            Collection<File> c = FileUtils.listFiles(directory, fileTypes.toArray(fileTypesArray), true);
            Iterator itFL = c.iterator();
            for(File f: c){
                logger.info(f.toString());
                try{
    
                    ArrayList<String> cmd = new ArrayList<>();
                    cmd.add(this.getUnzipSoftware());
                    cmd.add("-p " + this.getNoOfThreads());
                    //-d is for decompress
//                    cmd.add("-d");
                    cmd.add(f.toString());
                    

                    /*
                    pigz -p 4 -d /data/ngsdata/project1/sra_data.fastq.gz 
                    */      

                    cmdZip = this.cleanPath(StringUtils.join(cmd, " "));
                    logger.info("pigz zip command:\t" + cmdZip);

                    Runtime rt = Runtime.getRuntime();
                    Process proc = rt.exec(cmdZip);
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
                catch(IOException|InterruptedException ex){
                    logger.error("error executing pigz unzip command\n" + cmdZip);
                    logger.error(ex.toString());
                    throw new IOException(STEP_ID_STRING + ": error executing pigz unzip command" + cmdZip);
                }
             
            }            
        }
        
        
        
    }
    
    
    
            
    @Override
    public void verifyInputData() throws IOException{
        
        logger.info("verify input data");  
        
        if(new File(this.getUnzipSoftware()).exists() == false){
            logger.error("unzip software not found at location < " + this.getUnzipSoftware() +">");
            throw new IOException("unzip software not found at location < " + this.getUnzipSoftware() +">");
        }
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

        configData.put(ID_SOFTWARE, "/usr/local/pigz");
        configData.put(ID_THREADS, 4);
        configData.put(ID_FILE_TYPES, new ArrayList<>(Arrays.asList("fastq", "fasta", "sam")));

        
        return configData;
    }
    
    
    

    @Override
    public void verifyOutputData(){
        
    }

    @Override
    public NGSStepSubclass getStepSubclass(){
        return STEP_SUBCLASS;
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

    /**
     * @return the fileTypes
     */
    public ArrayList<String> getFileTypes() {
        return fileTypes;
    }

    /**
     * @param fileTypes the fileTypes to set
     */
    public void setFileTypes(ArrayList<String> fileTypes) {
        this.fileTypes = fileTypes;
    }

//    @Override
//    public void executeSingle() throws IOException {
//        logger.info(STEP_ID_STRING + ": execute step");        
//        
//        String projectRoot = getStepInputData().getProjectRoot() + FILESEPARATOR + stepInputData.getProjectID();
//        File directory = new File(projectRoot);
//        String cmdZip = "";
//        for(String fileType: fileTypes){
//            String fileTypesArray[] = new String[fileTypes.size()];
//            Collection<File> c = FileUtils.listFiles(directory, fileTypes.toArray(fileTypesArray), true);
//            Iterator itFL = c.iterator();
//            for(File f: c){
//                logger.info(f.toString());
//                try{
//    
//                    ArrayList<String> cmd = new ArrayList<>();
//                    cmd.add(this.getUnzipSoftware());
//                    cmd.add("-p " + this.getNoOfThreads());
//                     //-d is for decompress
////                    cmd.add("-d");
//                    cmd.add(f.toString());
//                    
//
//                    /*
//                    pigz -p 4 -d /data/ngsdata/project1/sra_data.fastq.gz 
//                    */      
//
//                    cmdZip = this.cleanPath(StringUtils.join(cmd, " "));
//                    logger.info("pigz zip command:\t" + cmdZip);
//
//                    Runtime rt = Runtime.getRuntime();
//                    Process proc = rt.exec(cmdZip);
//                    BufferedReader brStdin  = new BufferedReader(new InputStreamReader(proc.getInputStream()));
//                    BufferedReader brStdErr = new BufferedReader(new InputStreamReader(proc.getErrorStream()));
//
//                        String line = null;
//                        logger.info("<OUTPUT>");
//                        while ( (line = brStdin.readLine()) != null)
//                            logger.info(line);
//                        logger.info("</OUTPUT>");
//
//                        logger.info("<ERROR>");
//                        while ( (line = brStdErr.readLine()) != null)
//                            logger.info(line);
//                        logger.info("</ERROR>");
//
//
//                    int exitVal = proc.waitFor();            
//                    logger.info("Process exitValue: " + exitVal);   
//
//                    brStdin.close();
//                    brStdErr.close();
//                }
//                catch(IOException|InterruptedException ex){
//                    logger.error("error executing pigz unzip command\n" + cmdZip);
//                    logger.error(ex.toString());
//                    throw new IOException(STEP_ID_STRING + ": error executing pigz unzip command" + cmdZip);
//                }
//             
//            }            
//        }
//    }

//    @Override
//    public void verifyInputDataSingle() throws IOException {
//        logger.info("verify input data");  
//        
//        if(new File(this.getUnzipSoftware()).exists() == false){
//            logger.error("unzip software not found at location < " + this.getUnzipSoftware() +">");
//            throw new IOException("unzip software not found at location < " + this.getUnzipSoftware() +">");
//        }
//    }
}
