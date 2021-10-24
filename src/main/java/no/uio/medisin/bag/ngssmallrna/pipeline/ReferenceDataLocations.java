/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * contains information about root paths to the different data types
 * 
 * @author simonray
 */
public class ReferenceDataLocations {
    
    static Logger                       logger                     = LogManager.getLogger();

    public final static     String ID_CONFIG_ID                    = "ReferenceData";
    public final static     String ID_GENOME_FOLDER                = "genomeRootFolder";
    public final static     String ID_MIRBASE_FOLDER               = "mirbaseFolder";
    public final static     String ID_TARGETSCAN_FOLDER            = "targetscanFolder";
    public final static     String ID_MIRGENEDB_FOLDER             = "mirgenedbFolder";
    
    public final static     String ID_REL_BOWTIE_PATH              = "Sequence/BowtieIndex/genome";
    public final static     String ID_REL_BOWTIE2_PATH             = "Sequence/Bowtie2Index/genome";
    public final static     String ID_REL_HISAT_PATH               = "Sequence/HisatIndex/genome";
    public final static     String ID_REL_ABUN_DATA_PATH           = "Sequence/AbundantSequences/abundantIndex1";
    public final static     String ID_REL_ABUN_BOWTIE2_DATA_PATH   = "Sequence/AbundantSequences/abundantIndex2";
    public final static     String ID_REL_WHOLE_GENSEQ_PATH        = "Sequence/WholeGenomeFasta";
    public final static     String ID_REL_WHOLE_GENSEQ_FA          = "Sequence/WholeGenomeFasta/genome.fa";
    public final static     String ID_GENOME_FEATURELIST_FILE      = "feature.list";
    public final static     String ID_GENE_ANNOTATION              = "/Annotation/Genes";
    public final static     String ID_REL_RIBOSOMAL_RNA            = "Sequence/ribosomalRNA";
    
    public final static     String ID_TSCAN_MIRFAMILY_FILE         = "miR_Family_Info.txt";
    public final static     String ID_TSCAN_PREDICTIONS_FILE       = "Predicted_Targets_Info.txt";
    
    
    
    
    private                 String genomeRootFolder;
    private                 String mirbaseFolder;
    private                 String targetscanFolder;
    private                 String mirgenedbfolder;

    
    public ReferenceDataLocations(){
    
    }
    
    
    public ReferenceDataLocations(String genRootFolder, 
            String mirbRootFolder, 
            String tscanRootFolder, 
            String mirgRootFolder){
        genomeRootFolder    = genRootFolder;
        mirbaseFolder       = mirbRootFolder;
        targetscanFolder    = tscanRootFolder;
        mirgenedbfolder     = mirgRootFolder;
        
    }
       
    
    public ReferenceDataLocations(HashMap options){
        
        genomeRootFolder    = (String)options.get(ReferenceDataLocations.ID_GENOME_FOLDER);
        mirbaseFolder       = (String)options.get(ReferenceDataLocations.ID_MIRBASE_FOLDER);
        targetscanFolder    = (String)options.get(ReferenceDataLocations.ID_TARGETSCAN_FOLDER);
        mirgenedbfolder     = (String)options.get(ReferenceDataLocations.ID_MIRGENEDB_FOLDER);
    }
    
    
    /**
     * in this method we are simply checking that the configuration file 
     * has all the entries we need. We don't check if the values are acceptable,
     * that is the role of the @link{NGSStep} class.
     * 
     * @param configData
     * @throws Exception 
     */
    
    public void parseConfigurationData(HashMap configData) throws Exception{
        logger.info(this.getClass() + ": verify configuration data");
        
        if(configData.get(ID_GENOME_FOLDER)==null) {
            throw new NullPointerException("<" + ID_GENOME_FOLDER + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_MIRBASE_FOLDER)==null) {
            throw new NullPointerException("<" + ID_MIRBASE_FOLDER + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_TARGETSCAN_FOLDER)==null) {
            throw new NullPointerException("<" + ID_TARGETSCAN_FOLDER + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_MIRGENEDB_FOLDER)==null) {
            throw new NullPointerException("<" + ID_MIRGENEDB_FOLDER + "> : Missing Definition in Configuration File");
        }

        logger.info("passed");
    }
    
    
    
    /**
     * 
     * @throws IOException 
     */
    public void verifyReferenceData() throws IOException{
        if(new File(ID_GENOME_FOLDER).exists()==false){
            throw new IOException("root genome folder <" + ID_GENOME_FOLDER + "> not found");
        }
        if(new File(ID_MIRBASE_FOLDER).exists()==false){
            throw new IOException("root mirbase folder <" + ID_MIRBASE_FOLDER + "> not found");
        }
        if(new File(ID_TARGETSCAN_FOLDER).exists()==false){
            throw new IOException("root targetscan folder <" + ID_TARGETSCAN_FOLDER + "> not found");
        }
        if(new File(ID_MIRGENEDB_FOLDER).exists()==false){
            throw new IOException("root targetscan folder <" + ID_MIRGENEDB_FOLDER + "> not found");
        }
        
    }
    

    
    /**
     * generate sample configuration data
     * 
     * @return 
     */
    public HashMap generateExampleConfigurationData() {

        HashMap configData = new HashMap();
        configData.put(ID_GENOME_FOLDER, "/data/genomes");
        configData.put(ID_MIRBASE_FOLDER, "/data/mirbase");
        configData.put(ID_TARGETSCAN_FOLDER, "/data/targetscan");
        configData.put(ID_MIRGENEDB_FOLDER, "/data/targetscan");
    
        return configData;
    }
    
    
    
    
    /**
     * @return the genomeRootFolder
     */
    public String getGenomeRootFolder() {
        return genomeRootFolder;
    }

    /**
     * @param genomeRootFolder the genomeRootFolder to set
     */
    public void setGenomeRootFolder(String genomeRootFolder) {
        this.genomeRootFolder = genomeRootFolder;
    }

    /**
     * @return the mirbaseFolder
     */
    public String getMirbaseFolder() {
        return mirbaseFolder;
    }

    /**
     * @param mirbaseFolder the mirbaseFolder to set
     */
    public void setMirbaseFolder(String mirbaseFolder) {
        this.mirbaseFolder = mirbaseFolder;
    }

    /**
     * @return the targetbaseFolder
     */
    public String getTargetscanFolder() {
        return targetscanFolder;
    }

    /**
     * @param targetbaseFolder the targetbaseFolder to set
     */
    public void setTargetscanFolder(String targetbaseFolder) {
        this.targetscanFolder = targetbaseFolder;
    }

    /**
     * @return the mirgenedbfolder
     */
    public String getMirgenedbfolder() {
        return mirgenedbfolder;
    }

    /**
     * @param mirgenedbfolder the mirgenedbfolder to set
     */
    public void setMirgenedbfolder(String mirgenedbfolder) {
        this.mirgenedbfolder = mirgenedbfolder;
    }
    
}
