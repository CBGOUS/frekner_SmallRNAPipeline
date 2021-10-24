/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

import org.testng.Assert;
import org.testng.annotations.Test;

/**
 *
 * @author simonray
 */
public class PipelineDataNGTest {
    
    public PipelineDataNGTest() {
    }

    @Test
    public static void missingPipelineName(){
        PipelineData pipelineData = new PipelineData(null, "projectID", "projectRoot", "dataRoot");
        Assert.assertEquals(Boolean.FALSE,  pipelineData.isDataComplete());
    }

    @Test
    public static void missingProjectID(){
        PipelineData pipelineData = new PipelineData("projectName", null, "projectRoot", "dataRoot");
        Assert.assertEquals(Boolean.FALSE,  pipelineData.isDataComplete());
    }

    @Test
    public static void missingProjectRoot(){
        PipelineData pipelineData = new PipelineData("projectName", "projectID", null, "dataRoot");
        Assert.assertEquals(Boolean.FALSE,  pipelineData.isDataComplete());
    }

    @Test
    public static void missingDataRoot(){
        PipelineData pipelineData = new PipelineData("projectName", "projectID", "projectRoot", null);
        Assert.assertEquals(Boolean.FALSE,  pipelineData.isDataComplete());
    }

    
}
