/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import java.util.ArrayList;
import no.uio.medisin.bag.ngssmallrna.pipeline.ReferenceDataLocations;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;

import org.apache.logging.log4j.Logger;

/**
 *  Specifies the information needed to locate the input data associated with a step
 * ProcessBuilder
 * @author sr
 */
public class InputDataForStep {
    
    static  Logger  logger = LogManager.getLogger();
    static  String  FileSeparator = System.getProperty("file.separator");

    private String                      projectID;
    private String                      projectRoot;
    private ReferenceDataLocations      dataLocations;
    private String                      inputFolder;
    private String                      outputFolder;
    private String                      parameterString;
    private ArrayList<SampleDataEntry>  sampleData;
    
    
    public InputDataForStep(String pid, String pRoot, 
            ReferenceDataLocations dataLoc, String inFolder, String outFolder, String paramString,
            ArrayList<SampleDataEntry>sdata){
        
        projectID       = pid;
        projectRoot     = pRoot;
        dataLocations   = dataLoc;
        inputFolder     = inFolder;
        outputFolder    = outFolder;
        parameterString = paramString;
        sampleData      = sdata;
        
    }
    
    
    /**
     * check the input data is good
     * @throws ValidationException 
     */
    public void verifyInputData() throws ParseException
    {
        logger.info("verify input data");
        
        if(projectRoot==null || projectRoot.isEmpty())
            throw new ParseException("missing projectRoot definition");
        // check the project root exists
        // check the project folder exists within root
        
        
        /*
        if (!(new File(projectRoot)).exists()){
            throw new IOException();
        }
        logger.info("Project root <" + projectRoot + "> exists");
        
        String projectFolder = projectRoot + System.getProperty("file.separator") + projectID;
        projectFolder = projectFolder.replace(FileSeparator + FileSeparator, FileSeparator).trim();
        if (!(new File(projectFolder)).exists()){
            logger.error("project folder <" + projectFolder + "> not found");
            throw new IOException();
        }
        logger.info("Project folder <" + projectFolder + ">  exists");
        
        inputFolder = projectFolder + System.getProperty("file.separator") + inputFolder;
        if (!(new File(inputFolder)).exists()){
            logger.error("input folder <" + inputFolder + "> not found");
            throw new IOException();
        }
        */        

        
    }
    /**
     * @return the projectID
     */
    public String getProjectID() {
        return projectID;
    }

    /**
     * @param projectID the projectID to set
     */
    public void setProjectID(String projectID) {
        this.projectID = projectID;
    }

    /**
     * @return the projectRoot
     */
    public String getProjectRoot() {
        return projectRoot;
    }

    /**
     * @param projectRoot the projectRoot to set
     */
    public void setProjectRoot(String projectRoot) {
        this.projectRoot = projectRoot;
    }

    /**
     * @return the sampleData
     */
    public ArrayList<SampleDataEntry> getSampleData() {
        return sampleData;
    }

    /**
     * @return the inputFolder
     */
    public String getInputFolder() {
        return inputFolder;
    }

    /**
     * @param inputFolder the inputFolder to set
     */
    public void setInputFolder(String inputFolder) {
        this.inputFolder = inputFolder;
    }

    /**
     * @return the outputFolder
     */
    public String getOutputFolder() {
        return outputFolder;
    }

    /**
     * @return the dataLocations
     */
    public ReferenceDataLocations getDataLocations() {
        return dataLocations;
    }

  /**
   * @return the parameterString
   */
  public String getParameterString() {
    return parameterString;
  }

  /**
   * @param parameterString the parameterString to set
   */
  public void setParameterString(String parameterString) {
    this.parameterString = parameterString;
  }
    
    
}
