/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import java.io.IOException;
import java.util.HashMap;

/**
 *
 * @author simonray
 */
public interface NGSBase {
    void                parseStepParameters() throws Exception;
    void                verifyInputData() throws IOException;
    void                verifyOutputData() throws  IOException;
    void                parseConfigurationData(HashMap configData) throws Exception;
    HashMap             generateExampleConfigurationData();
    void                execute() throws IOException;
    String              longStepDescription();
    String              shortStepDescription();
    NGSStepSubclass     getStepSubclass();
   
}
