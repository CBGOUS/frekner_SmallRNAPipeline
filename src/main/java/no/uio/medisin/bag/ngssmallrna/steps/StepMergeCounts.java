/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.HashMap;

import no.uio.medisin.bag.core.mirna.MiRNA;
import no.uio.medisin.bag.core.mirna.MiRNASet;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


/**
 *   1. parse SAM file to extract and process the miRNA reads to determine isomiR content
 *   2. merge the counts from each sample to generate a single count file 
 *      that can be used for differential expression analysis
 * 
 *   Input is a SAM file
 * 
 * @author sr
 */

public class StepMergeCounts extends NGSStep{
    
    static Logger                       logger                      = LogManager.getLogger();
    
    public static final NGSStepSubclass STEP_SUBCLASS               = NGSStepSubclass.DATABLE;

    public  static final String     STEP_ID_STRING                  = "MergeCounts";
    private static final String     ID_REF_GENOME                   = "host";
    private static final String     ID_MIRBASE_VERSION              = "mirbaseVersion";

    private static final String     MIR_COUNTS_EXTENSION            = ".trim.clp.gen.mircounts.tsv";
    private static final String     RAW_INPUT_EXTENSION             = ".fastq.gz";
    private String                  miRBaseRelease                  = "20";
    private String                  ReferenceGenome                 = "";
    private double                  pValue                          = 0.05;
        
        
    private String                  mergedCountsFile;
    
    
    private MiRNASet                miRBaseMiRNAList                = new MiRNASet();


    public StepMergeCounts(){
        classSubtype = NGSStepSubclass.NOTDATABLE;
    }

    /**
     * 
     * @param sid 
     */
    public StepMergeCounts(InputDataForStep sid){
        classSubtype = NGSStepSubclass.NOTDATABLE;
        stepInputData = sid;
    }


    @Override
    public String shortStepDescription(){
      return "generate single merged count file from individual files";
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "generate single merged count file from individual files\n";
    }


    @Override
    public void parseStepParameters() throws Exception{
      
    }
    


    
    /**
     * in this method we are simply checking that the configuration file 
     * has all the entries we need. We dont check if the values are acceptable
     * that is the role of the NGSStep.
     * 
     * @param configData
     * @throws Exception 
     */
    @Override
    public void parseConfigurationData(HashMap configData) throws Exception{
        
        logger.info(STEP_ID_STRING + ": verify configuration data");
        
        if(configData.get(ID_REF_GENOME)==null) {
            logger.error("<" + ID_REF_GENOME + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_REF_GENOME + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_MIRBASE_VERSION)==null) {
            logger.error("<" + ID_MIRBASE_VERSION + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_MIRBASE_VERSION + "> : Missing Definition in Configuration File");
        }

        /*try{
            this.setMiRBaseRelease((Integer) configData.get(ID_MIRBASE_VERSION));
        }
        catch(Exception exNm){
            logger.error(ID_MIRBASE_VERSION + " <" + configData.get(ID_MIRBASE_VERSION) + "> is not an integer");
            throw new NumberFormatException(ID_MIRBASE_VERSION + " <" + configData.get(ID_MIRBASE_VERSION) + "> is not an integer");
        }        
        if (this.getMiRBaseRelease() <= 0){
            logger.error(ID_MIRBASE_VERSION + " <" + configData.get(ID_MIRBASE_VERSION) + "> must be positive integer");
            throw new IllegalArgumentException(ID_MIRBASE_VERSION + " <" + configData.get(ID_MIRBASE_VERSION) + "> must be positive integer");
        }
        */
        this.setMiRBaseRelease((String)configData.get(ID_MIRBASE_VERSION));
       
        this.setReferenceGenome((String) configData.get(ID_REF_GENOME));
        if(this.getReferenceGenome().length() !=3 ){
            logger.error(ID_REF_GENOME + " <" + configData.get(ID_REF_GENOME) + "> must be a 3 letter string");            
            throw new IllegalArgumentException(ID_REF_GENOME + " <" + configData.get(ID_REF_GENOME) + "> must be a 3 letter string");            
        }

        logger.info("passed");
    }
    
    
    
    
    
    
    @Override
    public void execute() throws IOException{
        
        logger.info(STEP_ID_STRING + ": execute step");        
        
        /*
            1. read in all sample count files and merge
            2. output merged count file
        */
        Boolean fA = new File(outFolder).mkdir();       
        if (fA) logger.info("created output folder <" + outFolder + "> for results" );
        
        String gffFileMirBase = this.cleanPath(getStepInputData().getDataLocations().getMirbaseFolder() 
                + FILESEPARATOR + this.getMiRBaseRelease() + FILESEPARATOR + this.getReferenceGenome() + ".gff3");
        String faFileMirBase = gffFileMirBase.replace("gff3", "mature.fa");
        miRBaseMiRNAList.loadMiRBaseData(this.getReferenceGenome(), gffFileMirBase, faFileMirBase);

        logger.info("Merging Count Files");
        String headerLine = "name";
        String[] countStrings = new String[miRBaseMiRNAList.getNumberOfEntries()];
        Arrays.fill(countStrings, "");
        
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            headerLine = headerLine.concat("\t" + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, ""));
            String  miRCountsFile  = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, MIR_COUNTS_EXTENSION));
            try{
                int m=0;
                BufferedReader brmiRCounts  = new BufferedReader(new FileReader(new File(miRCountsFile)));
                    String countLine = "";
                    while((countLine = brmiRCounts.readLine()) != null){
                        logger.info(Integer.toString(m )+ ":" + countLine);
                        countStrings[m] = countStrings[m].concat("\t" + countLine.split("\t")[1].trim() );
                        m++;
                    }
                brmiRCounts.close();
            }
            catch(Exception ex){
                logger.error("error reading count files for merging <" + miRCountsFile + "> \n" + ex.toString());
                throw new IOException(STEP_ID_STRING + " :error reading counts file  <" + miRCountsFile + ">");
            }
        }
        
        logger.info("Writing merged count files");
        mergedCountsFile        = outFolder + FILESEPARATOR + getStepInputData().getProjectID() + ".merged.mirna_counts.tsv"; 
        
        try{
            BufferedWriter bwMC = new BufferedWriter(new FileWriter(new File(mergedCountsFile)));
            
            bwMC.write(headerLine + "\n");
            int m=0;
            for(MiRNA miR: miRBaseMiRNAList.getMiRBaseMiRNAList()){
                //bwMC.write(miR.getMimatID() + "|" + miR.getMimatAlias()+ countStrings[m] + "\n");
                bwMC.write(miR.getMimatID() + "|" + miR.getName()+ countStrings[m] + "\n");
                m++;
            }
            
            bwMC.close();
        }
        catch(IOException exIO){
            logger.info("error writing merged counts File <" + mergedCountsFile + ">\n" + exIO);        
            throw new IOException(STEP_ID_STRING + " :error writing merged counts file  <" + mergedCountsFile + ">");
        }
        
        
        logger.info(STEP_ID_STRING + ": completed");
        
    }
    
    
    
    /**
     * Verify Input Data 
     * 
     * @throws IOException
     */        
    @Override
    public void verifyInputData() throws IOException{
        
        logger.info("verify input data");        
        this.setPaths();
        
//        String gffFileMirBase = this.cleanPath(getStepInputData().getDataLocations().getMirbaseFolder() 
//                + FILESEPARATOR + this.getMiRBaseRelease() + FILESEPARATOR + this.getReferenceGenome() + ".gff3");
//        if (new File(gffFileMirBase).exists()==false){
//            logger.error("no annotation file was found for mirBase HOST:<" 
//                    + this.getReferenceGenome() + "> VERSION: <"+ this.getMiRBaseRelease() + "> at location <" 
//                    + gffFileMirBase + ">");
//            throw new IOException("no annotation file was found for mirBase HOST:<" 
//                    + this.getReferenceGenome() + "> VERSION: <"+ this.getMiRBaseRelease() + "> at location <" 
//                    + gffFileMirBase + ">");
//        }
//                
//        String faFileMirBase = gffFileMirBase.replace("gff3", "mature.fa");
//        if (new File(gffFileMirBase).exists()==false){
//            logger.error("no fasta file was found for mirBase HOST:<" 
//                    + this.getReferenceGenome() + "> VERSION: <"+ this.getMiRBaseRelease() + "> at location <" 
//                    + faFileMirBase + ">");
//            throw new IOException("no fasta file was found for mirBase HOST:<" 
//                    + this.getReferenceGenome() + "> VERSION: <"+ this.getMiRBaseRelease() + "> at location <" 
//                    + faFileMirBase + ">");
//        }
                
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            String  miRCountsFile  = inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, MIR_COUNTS_EXTENSION);            
            if ((new File(miRCountsFile)).exists()==false){
                throw new IOException(STEP_ID_STRING + ": counts file <" + sampleData + "> does not exist");
            }
            
        }
    }
    
    
    
    
    /**
     * generate sample configuration data so the user can see what can be
     * specified
     *
     * @return
     */
    @Override
    public HashMap generateExampleConfigurationData() {

        logger.info(STEP_ID_STRING + ": generate example configuration data");

        HashMap configData = new HashMap();

        configData.put(ID_REF_GENOME, "hsa");
        configData.put(ID_MIRBASE_VERSION, 20);

       
        return configData;
    }


    
    
    
    @Override
    public void verifyOutputData(){
        
    }


    @Override
    public NGSStepSubclass getStepSubclass(){
        return STEP_SUBCLASS;
    }
    
    


    /**
     * @return the pValue
     */
    public double getpValue() {
        return pValue;
    }

    /**
     * @param pValue the pValue to set
     */
    public void setpValue(double pValue) {
        this.pValue = pValue;
    }

    /**
     * @return the miRBaseRelease
     */
    public String getMiRBaseRelease() {
        return miRBaseRelease;
    }

    /**
     * @param miRBaseRelease the miRBaseRelease to set
     */
    public void setMiRBaseRelease(String miRBaseRelease) {
        this.miRBaseRelease = miRBaseRelease;
    }

    /**
     * @return the ReferenceGenome
     */
    public String getReferenceGenome() {
        return ReferenceGenome;
    }

    /**
     * @param ReferenceGenome the ReferenceGenome to set
     */
    public void setReferenceGenome(String ReferenceGenome) {
        this.ReferenceGenome = ReferenceGenome;
    }


}
