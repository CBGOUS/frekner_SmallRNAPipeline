/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;


import org.apache.logging.log4j.Logger;


/**
 * THIS CLASS SHOULD NOT BE CHANGED
 * @author sr
 */
abstract public class NGSStep implements NGSBase{
    
    protected               InputDataForStep    stepInputData   = null;
    protected               ResultDataForStep   stepResultData  = null;
    
    protected               String              inFolder        = null;
    protected               String              outFolder       = null;
    
    protected static final  String              FILESEPARATOR   = System.getProperty("file.separator");
    protected static final  String              DISCOUNT_PARAMETER_LIMIT     = "NA";
    
    protected static final  String              BOOLEAN_NAME    = Boolean.class.getName();
    public                  NGSStepSubclass     classSubtype;
    
    public static final     int                 LOGINDENT       = 6;
    
    

   

/*    
    abstract void       verifyInputData() throws IOException, NullPointerException;
    abstract void       verifyOutputData();
    abstract void       parseConfigurationData(HashMap configData) throws Exception;
    abstract HashMap    generateExampleConfigurationData();
    abstract void       execute() throws IOException;
*/    
    
    /**
     * set paths for input and output data folders
     * 
     */
    final void setPaths(){
        
        //getStepInputData().getInputFolder();
        String projectFolder = getStepInputData().getProjectRoot() + FILESEPARATOR + getStepInputData().getProjectID();
        projectFolder = projectFolder.replace(FILESEPARATOR + FILESEPARATOR, FILESEPARATOR).trim();
        
        inFolder = projectFolder + FILESEPARATOR + getStepInputData().getInputFolder();
        outFolder = projectFolder + FILESEPARATOR + getStepInputData().getOutputFolder();
        
    }
    
    /**
     * strip out duplicate folder delimiter from a file path
     * 
     * @param path
     * @return 
     */
    final String cleanPath(String path){
        return path.replace(FILESEPARATOR + FILESEPARATOR, FILESEPARATOR);
    }
    
    
    
    
 
    /**
     * Check the specified parameter is valid
     * For numerical parameters this requires not null and within the specified 
     * upper and lower limits
     * 
     * @param className
     * @param paramID
     * @param param
     * @param lowerLimit
     * @param upperLimit
     * @param logger
     * @return
     * @throws Exception 
     */
    protected String checkParameter(String className, String paramID, String param, String lowerLimit, String upperLimit, Logger logger) throws Exception{
        int             iVal;
        double          dVal;
        Boolean         bVal;
        String          sVal;
                
        if(Integer.class.getName().contains(className)){
            iVal = checkIntegerParameter(paramID, param, lowerLimit, upperLimit, logger);
            return Integer.toString(iVal);
        }
        
        if(Double.class.getName().contains(className)){
            dVal = checkDoubleParameter(paramID, param, lowerLimit, upperLimit, logger);
            return Double.toString(dVal);
        }
        
        if(Boolean.class.getName().contains(className)){
            bVal = checkBooleanParameter(paramID, param, lowerLimit, upperLimit, logger);
            return Boolean.toString(bVal);
        }
        
        if(String.class.getName().contains(className)){
            sVal = this.checkStringParameter(paramID, param, lowerLimit, upperLimit, logger);
            return sVal;
        }
                
        return null;
    }
    
    
    
    
    
    /**
     * 
     * check the specified Integer parameter is valid, i.e., not null and within the specified range
     * if there is only an upper or lower limit, the other limit can be set to NGSStep.DISCOUNT_PARAMETER_LIMIT
     * 
     * @param paramID
     * @param param
     * @param lowerLimit
     * @param upperLimit
     * @param logger
     * @return 
     */
    private int checkIntegerParameter(String paramID, String param, String lowerLimit, String upperLimit, Logger logger) throws Exception{
        
        
        int             iVal;
        
        
        try{
            iVal = Integer.parseInt(param);
        }
        catch(NumberFormatException exNm){
            logger.info(paramID + " <" + param + "> is not an integer");
            logger.error(paramID + " <" + param + "> is not an integer");
            throw new NumberFormatException(paramID + " <" + param + "> is not an integer");
        }        

        if(lowerLimit.equals(NGSStep.DISCOUNT_PARAMETER_LIMIT) && upperLimit.equals(NGSStep.DISCOUNT_PARAMETER_LIMIT))
            return iVal;

        if(upperLimit.equals(NGSStep.DISCOUNT_PARAMETER_LIMIT)){
            if (iVal < Integer.parseInt(lowerLimit)){
                logger.info(paramID + " <" + param + "> must be >= " + lowerLimit);
                logger.error(paramID + " <" + param + "> must be >= " + lowerLimit);
                throw new IllegalArgumentException(paramID + " <" + param + "> must be >= " + lowerLimit);
            }
            else 
                return iVal;
        }

        if(lowerLimit.equals(NGSStep.DISCOUNT_PARAMETER_LIMIT)){
            if (iVal > Integer.parseInt(upperLimit)){
                logger.info(paramID + " <" + param + "> must be > " + upperLimit);
                logger.error(paramID + " <" + param + "> must be > " + upperLimit);
                throw new IllegalArgumentException(paramID+ " <" + param + "> must be > " + upperLimit);
            }
            else 
                return iVal;
        }

        if(iVal <Integer.parseInt(lowerLimit) 
                || iVal > Integer.parseInt(upperLimit)){
                logger.info(paramID + " <" + param 
                        + "> must be must be >= " + lowerLimit + " <= " + upperLimit);
                logger.error(paramID + " <" + param 
                        + "> must be must be >= " + lowerLimit + " <= " + upperLimit);
                throw new IllegalArgumentException(paramID + " <" + param 
                        + "> must be must be >= " + lowerLimit + " <= " + upperLimit);

        }
        else{
            return iVal;            
        }
        
    }
    
    
    
    
    
    
    /**
     * 
     * check the specified Double parameter is valid, i.e., not null and within the specified range
     * if there is only an upper or lower limit, the other limit can be set to NGSStep.DISCOUNT_PARAMETER_LIMIT
     * 
     * @param paramID
     * @param param
     * @param lowerLimit
     * @param upperLimit
     * @param logger
     * @return 
     */
    private double checkDoubleParameter(String paramID, String param, String lowerLimit, String upperLimit, Logger logger) throws Exception{
        
        
        double             dVal;
        try{
            dVal = Double.parseDouble(param);
        }
        catch(NumberFormatException exNm){
            logger.info(paramID + " <" + param + "> is not an integer");
            logger.error(paramID + " <" + param + "> is not an integer");
            throw new NumberFormatException(paramID + " <" + param + "> is not an integer");
        }        

        if(lowerLimit.equals(NGSStep.DISCOUNT_PARAMETER_LIMIT) && upperLimit.equals(NGSStep.DISCOUNT_PARAMETER_LIMIT))
            return dVal;

        if(upperLimit.equals(NGSStep.DISCOUNT_PARAMETER_LIMIT)){
            if (dVal < Double.parseDouble(lowerLimit)){
                logger.info(paramID + " <" + param + "> must be >= " + lowerLimit);
                logger.error(paramID + " <" + param + "> must be >= " + lowerLimit);
                throw new IllegalArgumentException(paramID + " <" + param + "> must be >= " + lowerLimit);
            }
            else 
                return dVal;
        }

        if(lowerLimit.equals(NGSStep.DISCOUNT_PARAMETER_LIMIT)){
            if (dVal > Double.parseDouble(upperLimit)){
                logger.info(paramID + " <" + param + "> must be > " + upperLimit);
                logger.error(paramID + " <" + param + "> must be > " + upperLimit);
                throw new IllegalArgumentException(paramID+ " <" + param + "> must be > " + upperLimit);
            }
            else 
                return dVal;
        }

        if(dVal <Double.parseDouble(lowerLimit) 
                || dVal > Double.parseDouble(upperLimit)){
                logger.info(paramID + " <" + param 
                        + "> must be must be >= " + lowerLimit + " <= " + upperLimit);
                logger.error(paramID + " <" + param 
                        + "> must be must be >= " + lowerLimit + " <= " + upperLimit);
                throw new IllegalArgumentException(paramID + " <" + param + "> must be > " + upperLimit);

        }
        else{
            return dVal;            
        }
        
    }
    
    
    
    
    
    
    /**
     * check the specified Boolean parameter is valid. i.e., can it be cast as a Boolean
     * 
     * @param paramID
     * @param param
     * @param lowerLimit
     * @param upperLimit
     * @param logger
     * @return
     * @throws Exception 
     */
    private boolean checkBooleanParameter(String paramID, String param, String lowerLimit, String upperLimit, Logger logger) throws Exception{
        
        
        boolean             bVal;
        try{
            bVal = Boolean.parseBoolean(param);
        }
        catch(NumberFormatException exNm){
            logger.info(paramID + " <" + param + "> cannot be cast as Boolean");
            logger.error(paramID + " <" + param + "> cannot be cast as Boolean");
            throw new NumberFormatException(paramID + " <" + param + "> cannot be cast as Boolean");
        }        
            
        return bVal;            
            
    }

    
    /**
     * check the specified String parameter is valid, i.e., not null and within the specified length range
     * if there is only an upper or lower limit, the other limit can be set to NGSStep.DISCOUNT_PARAMETER_LIMIT
     * 
     * @param paramID
     * @param param
     * @param lowerLimit
     * @param upperLimit
     * @param logger
     * @return
     * @throws Exception 
     */
    private String checkStringParameter(String paramID, String param, String lowerLimit, String upperLimit, Logger logger) throws Exception{
        
        
        int             iVal;
        
        
        try{
            iVal = Integer.parseInt(param);
        }
        catch(NumberFormatException exNm){
            logger.info(paramID + " <" + param + "> is not an integer");
            logger.error(paramID + " <" + param + "> is not an integer");
            throw new NumberFormatException(paramID + " <" + param + "> is not an integer");
        }        

        if(lowerLimit.equals(NGSStep.DISCOUNT_PARAMETER_LIMIT) && upperLimit.equals(NGSStep.DISCOUNT_PARAMETER_LIMIT))
            return param;

        if(upperLimit.equals(NGSStep.DISCOUNT_PARAMETER_LIMIT)){
            if (iVal < Integer.parseInt(lowerLimit)){
                logger.info(paramID + " <" + param + "> must be >= " + lowerLimit);
                logger.error(paramID + " <" + param + "> must be >= " + lowerLimit);
                throw new IllegalArgumentException(paramID + " <" + param + "> must be >= " + lowerLimit);
            }
            else 
                return param;
        }

        if(lowerLimit.equals(NGSStep.DISCOUNT_PARAMETER_LIMIT)){
            if (iVal > Integer.parseInt(upperLimit)){
                logger.info(paramID + " <" + param + "> must be > " + upperLimit);
                logger.error(paramID + " <" + param + "> must be > " + upperLimit);
                throw new IllegalArgumentException(paramID+ " <" + param + "> must be > " + upperLimit);
            }
            else 
                return param;
        }

        if(iVal <Integer.parseInt(lowerLimit) 
                || iVal > Integer.parseInt(upperLimit)){
                logger.info(paramID + " <" + param 
                        + "> must be must be >= " + lowerLimit + " <= " + upperLimit);
                logger.error(paramID + " <" + param 
                        + "> must be must be >= " + lowerLimit + " <= " + upperLimit);
                throw new IllegalArgumentException(paramID + " <" + param 
                        + "> must be must be >= " + lowerLimit + " <= " + upperLimit);

        }
        else{
            return param;            
        }
        
            
    }

    /**
     * @return the stepInputData
     */
    public InputDataForStep getStepInputData() {
        return stepInputData;
    }
    
    
    public void setStepInputData(InputDataForStep sid){
        stepInputData = sid; 
    }
}
