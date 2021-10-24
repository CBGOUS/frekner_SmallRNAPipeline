/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

import java.io.IOException;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author sr
 */
public class SmallNGSCmd {
    
    static Logger logger = LogManager.getRootLogger();

    static Options options = new Options();
    
    public static void main(String args[]) throws Exception{
        
        SmallNGSPipeline smallPipeline = new SmallNGSPipeline();
        parseArguments(args, smallPipeline);
        
        if(smallPipeline.getGenerateExampleConfigurationFile()){
            try{smallPipeline.generateSampleConfigurationFile();}
            catch(IOException exIO){
                logger.error(exIO.toString());
            }
            return;

        }
        try{
            smallPipeline.preparePipeline();
            smallPipeline.prepareSteps();
            smallPipeline.executePipeline();
        }
        catch(IOException exIO){
            logger.info("IOException while executing pipeline");
            logger.info(exIO.toString());
            logger.error(exIO);
        }
        catch(Exception exEx){
            logger.info("Exception while executing pipeline");
            logger.info(exEx.toString());
            logger.error(exEx);
        }
        
        
    }
    
    
    
    /**
     * 
     * get pipeline settings
     * User must provide the following
     *   run config file    : all the information needed to run the pipeline (e.g. data folders)
     *   data file          : contains information about samples and grouping
     *   steps file         : contains information about what type of analyses should be performed
     * 
     * @param args 
     * @param smallPipeline 
     */
    public static void parseArguments(String args[], SmallNGSPipeline smallPipeline) throws ParseException{
        
        
        logger.info("parse arguments");
        options.addOption("h", "help",             false,  "view help");
        options.addOption("c", "config",           true,   "run configuration file");
        options.addOption("p", "pipeline",         true,   "pipeline configuration file");
        options.addOption("d", "data",             true,   "data file list with grouping");
        options.addOption("G", "generate",         true,   "generate sample run files");
        options.addOption("D", "execute_by_data",  false,  "execute by data");
        options.addOption("S", "execute_by_steps", false,  "execute by steps");
        
        CommandLineParser clParser = new BasicParser();
        CommandLine cmd = null;
        
        
        try{

          
            cmd = clParser.parse(options, args);
                
            if(cmd.hasOption("h")){
                printHelp();
            }
            
            if (cmd.hasOption("G")) {
                logger.info("generate example run files ");
                smallPipeline.setConfigurationFile(cmd.getOptionValue("G"));
                smallPipeline.setGenerateExampleConfigurationFile(Boolean.TRUE);
            }   
            
            if (cmd.hasOption("D")){
                logger.info("execute program by data set");
                smallPipeline.setExecuteModel(SmallNGSPipeline.EXECUTE_BY_DATA);
            }
            
            if (cmd.hasOption("S")){
                logger.info("execute program by steps");
                smallPipeline.setExecuteModel(SmallNGSPipeline.EXECUTE_BY_STEPS);
            }
            
            if (cmd.hasOption("c")) {
                logger.info("run configuration file set to " + cmd.getOptionValue("c"));
                smallPipeline.setConfigurationFile(cmd.getOptionValue("c"));
            }                    
            else
                throw new ParseException("no configuration file was specified") ;
            
            if (cmd.hasOption("p")) {
                logger.info("pipeline file set to " + cmd.getOptionValue("p"));
                smallPipeline.setPipelineFile(cmd.getOptionValue("p"));
            }                    
            else
                throw new ParseException("no pipeline file was specified") ;
            
            if (cmd.hasOption("d")) {
                logger.info("NGS data filelist set to " + cmd.getOptionValue("d"));
                smallPipeline.setDataFile(cmd.getOptionValue("d"));
            }                    
            else
                throw new ParseException("no NGS data filelist was specified") ;

        }
        catch(ParseException exPa){
            logger.info("error parsing parameters:" + System.lineSeparator() + exPa);
            logger.error("error parsing parameters::" + System.lineSeparator() + exPa);
            throw new ParseException("error parsing parameters");
        }
    }
    

    
    
    /**
     * print command line usage
     * 
     */
    public static void printHelp(){
        printBanner();
        HelpFormatter formatter = new HelpFormatter();

        formatter.printHelp("command line options", options);
        
        System.exit(0);
        
    }
    
    
    
    /**
     * print program info
     * 
     */
    public static void printBanner(){
        logger.info("\n\n\n"
                + "    =======================================================================\n"
                + "    |  NGS smallRNA pipeline:                                             |\n"
                + "    |    performs a range of analyses on smallRNA NGS data                |\n"
                + "    =======================================================================\n\n\n");
        
        logger.info("*** report bugs to simon.rayner@medisin.uio.no\n");
    }
    
    
}
