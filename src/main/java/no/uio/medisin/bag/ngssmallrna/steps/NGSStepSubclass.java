/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

/**
 * specifies a key property of a NGSStep. 
 * 
 * DATABLE means the step can be subdivided
 * into substeps, one for each sample to produce a smaller hard disk footprint.
 * 
 * NOTDATABLE means the step has remain intact and all samples have to be processed
 * within a single instance of the step (e.g., for differential expression analysis,
 * all input data from all samples must be available for the analysis to proceed.
 * 
 * @author simonray
 */
public enum NGSStepSubclass {
    DATABLE, NOTDATABLE
}
