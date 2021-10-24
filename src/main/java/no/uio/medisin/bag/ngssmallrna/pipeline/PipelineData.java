/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

import java.util.ArrayList;
import java.util.List;
import no.uio.medisin.bag.ngssmallrna.steps.NGSRunStepJSONData;

/**
 * contains the project level data for the pipeline.
 * The individual steps need this information to access the data
 * 
 * @author sr
 */
public class PipelineData {
 
        private String pipelineName;
        private String projectID;
        private String projectRoot;
        private String dataRoot;
        
        private final List <NGSRunStepJSONData> stepsData = new ArrayList<>();
        
        public PipelineData(){
            
        }
        
        public PipelineData(String pName, String pID, String pRoot, String dRoot){
            pipelineName = pName;
            projectID = pID;
            projectRoot = pRoot;
            dataRoot = dRoot;
        }
 
	@Override
	public String toString() {
            String str = "User [" + pipelineName + "]\n";
                for(NGSRunStepJSONData step: getStepsJSONData()){
                    str += step.toString() + "\n";
                }
		return str;
	}

        
        
        
        /**
         * check that all values have been set - if any are missing then
         * complete paths to data cannot be formed
         * 
         * @return 
         */
        public Boolean isDataComplete() {
            return this.pipelineName != null
                    && this.projectID != null
                    && this.projectRoot != null
                    && this.dataRoot != null;
        }   
        
        
    /**
     * @return the steps
     */
    public List <NGSRunStepJSONData> getStepsJSONData() {
        return stepsData;
    }

    /**
     * @return the projectID
     */
    public String getProjectID() {
        return projectID;
    }

    /**
     * @return the projectRoot
     */
    public String getProjectRoot() {
        return projectRoot;
    }

    /**
     * @return the dataRoot
     */
    public String getDataRoot() {
        return dataRoot;
    }
}
