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

public class StepExit extends NGSStep {
    
    static Logger                       logger = LogManager.getLogger();
    
    public static final NGSStepSubclass STEP_SUBCLASS           = NGSStepSubclass.DATABLE;

    public static final String          STEP_ID_STRING          = "Exit";

    

    public StepExit(){
        classSubtype = NGSStepSubclass.DATABLE;
    }
    
    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepExit(InputDataForStep sid){
        classSubtype = NGSStepSubclass.DATABLE;
        stepInputData = sid;
    }
    
    
    @Override
    public void parseStepParameters() throws Exception{
      
    }
    

    @Override
    public String shortStepDescription(){
      return "Exit the pipeline";
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "Exit the pipeline.\n";
    }
    
    
    

    
    /**
     * This parses out the run parameters for this step
     * For this step. there won't be any data to validate
     * 
     * @param configData
     * @throws ValidationException 
     */
    @Override
    public void parseConfigurationData(HashMap configData) throws ParseException{

        logger.info(STEP_ID_STRING + ": verify configuration data");
       
 
        logger.info("passed");
    }
    
    
    
    /**
     * Nothing to execute for this step
     * @throws IOException 
     */
    @Override
    public void execute() throws IOException{
        /*
            sample fastqcSoftware command       
            fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam]  
         */
        logger.info(STEP_ID_STRING + ": execute step");


        logger.info(STEP_ID_STRING + ": completed");
    }
    
    
    
            
    /**
     * this should be called prior to executing the step.
     * no verification is necessary for this step
     * 
     * @throws IOException
     */
    @Override
    public void verifyInputData() throws IOException{
        
        logger.info("verify input data");        

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
     * There is no configuration data for the Exit step 
     * 
     * @return 
     */
    @Override
    public HashMap generateExampleConfigurationData() {

        logger.info(STEP_ID_STRING + ": generate example configuration data");
        
        HashMap<String, Object> configData = new HashMap();
        

        
        return configData;
        
    }

        
}
