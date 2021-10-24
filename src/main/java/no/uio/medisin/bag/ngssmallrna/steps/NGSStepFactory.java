/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import java.io.IOException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * StepFactory for creating new NGSSteps
 * 
 * @author simonrayner
 */
public class NGSStepFactory {
  static Logger logger = LogManager.getLogger();
  
  public static NGSStep getNewStep(String stepType) throws IOException{
    switch (stepType){
        
      case StepGetReadsIntersectingFeatures.STEP_ID_STRING:
        return new StepGetReadsIntersectingFeatures();

      case StepSingleReadAdapterTrim.STEP_ID_STRING:
        return new StepSingleReadAdapterTrim();

      case StepCollapseReads.STEP_ID_STRING:
        return new StepCollapseReads();

      case StepBowtieMapSingleReads.STEP_ID_STRING:
        return new StepBowtieMapSingleReads();
        
      case StepUnzipInputFiles.STEP_ID_STRING:
        return new StepUnzipInputFiles();
        
      case StepPairedEndTrimAdapters.STEP_ID_STRING:
        return new StepPairedEndTrimAdapters();

      case StepFastqQC.STEP_ID_STRING:
        return new StepFastqQC();
        
      case StepParseSAMForMiRNAs.STEP_ID_STRING:
        return new StepParseSAMForMiRNAs();
        
      case StepMergeCounts.STEP_ID_STRING:
        return new StepMergeCounts();
        
      case StepDEwithEdgeR.STEP_ID_STRING:
          return new StepDEwithEdgeR();
        
      case StepMergeIsomiRCounts.STEP_ID_STRING:
          return new StepMergeIsomiRCounts();
        
      case "exit":
          return null;

      default:
          logger.error("unknown step:");
          throw new IOException("unknown step: " + stepType);
          
    }
  }
}
