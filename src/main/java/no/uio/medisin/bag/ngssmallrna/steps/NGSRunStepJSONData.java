/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import org.apache.commons.cli.ParseException;


/**
 * This class stores the step data specified in the JSON project file
 * For example:
 * 
 *      {"step":"AnalyzeStartPositions", 
 *       "inputFolder":"inFld1", 
 *       "outputFolder":"outFold1", 
 *       "parameters": "--bed_file=bedfile1.bed, --bleed=3"},
 *
 * @author sr
 */
public class NGSRunStepJSONData {
    
    private String step;
    private String inputFolder;
    private String outputFolder;
    private String parameterString;

    
    /**
      this empty constructor is necessary for JSON serialization of data
    */
    public NGSRunStepJSONData(){
        
    }
    
    public NGSRunStepJSONData(String sType, String inFolder, String outFolder){
        step = sType;
        inputFolder = inFolder;
        outputFolder = outFolder;
    }
    
    
    /**
     * verify the specified data is complete
     * 
     * @throws ValidationException 
     */
    public void verifyData() throws ParseException{
        if(step == null || step.isEmpty()){
            throw new ParseException("missing Step name");
        }
        if(inputFolder == null || inputFolder.isEmpty()){
            throw new ParseException("missing InputFolder definition");
        }
        if(outputFolder == null || outputFolder.isEmpty()){
            throw new ParseException("missing OutputFolder definition");
        }
    }
    
    
    
    
    @Override
    public String toString() {
        return "step=" + step + "\t inputFolder:" + inputFolder + "\t outputFolder:" + outputFolder + "\t parameterString:" + parameterString;
    }    
    
    /**
     * @return the name
     */
    public String getStep() {
        return step;
    }

    /**
     * @param name the name to set
     */
    public void setStep(String name) {
        this.step = name;
    }

    /**
     * @return the inputFileList
     */
    public String getInputFolder() {
        return inputFolder;
    }

    /**
     * @param inputFileList the inputFileList to set
     */
    public void setInputFolder(String inputFileList) {
        this.inputFolder = inputFileList;
    }

    /**
     * @return the outputFileList
     */
    public String getOutputFolder() {
        return outputFolder;
    }

    /**
     * @param outputFileList the outputFileList to set
     */
    public void setOutputFolder(String outputFileList) {
        this.outputFolder = outputFileList;
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
