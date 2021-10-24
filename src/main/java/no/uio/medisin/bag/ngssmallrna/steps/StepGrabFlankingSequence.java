/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import no.uio.medisin.bag.core.reference.GFFEntry;
import no.uio.medisin.bag.core.reference.GFFSet;
import no.uio.medisin.bag.core.sequence.GenomeSeq;
import no.uio.medisin.bag.ngssmallrna.pipeline.ReferenceDataLocations;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.logging.log4j.LogManager;

import org.apache.logging.log4j.Logger;



/**
 *  Grab Flanking Sequence
 *  1. Input a list of sequences in FASTA format.
 *  2. Read the location information that is specified in the header line
 *  3. add the flanking sequence on either side
 *  4. output the new sequence with updated position information
 * 
 *  Header information is in the format:
 *    >msy-6088:5|chr4:273774-273790 6088:5|chr4:273774 msy 6088:5|chr4:273774
 *    GGTGGAATAGTACCTGT
 *  i.e. Name|Chromosome:Start or Name|Chromosome:Start-Stop
 * 
 * @author sr
 * 
 */

public class StepGrabFlankingSequence extends NGSStep{
    
    static  Logger                      logger = LogManager.getLogger();

    public static final NGSStepSubclass STEP_SUBCLASS               = NGSStepSubclass.DATABLE;

    public static final String          STEP_ID_STRING              = "GrabFlankingSequence";
    private static final String         ID_REF_GENOME               = "host";
    private static final String         FIVE_PRIME_FLANK_LEN        = "FivePrimeFlankLength";
    private static final String         THREE_PRIME_FLANK_LEN       = "ThreePrimeFlankLength";
    
    
    private static final String         INFILE_EXTENSION            = ".features.tsv";
    private static final String         OUTFILE_EXTENSION           = ".flanked.fasta";
    private static final String         MIR_COUNTS_EXTENSION        = ".trim.clp.gen.mircounts.tsv";
    
    private             String          rootDataFolder              = "";
    private             String          ReferenceGenome             = "";
    private             GenomeSeq       genomeFasta;
    private             GFFSet          featureSet;

    private Integer                     fiveFlankLen               = 0;
    private Integer                     threeFlankLen              = 0;
    
    
    

    public StepGrabFlankingSequence(){
        classSubtype = NGSStepSubclass.DATABLE;
    }
    
    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepGrabFlankingSequence(InputDataForStep sid){
        classSubtype = NGSStepSubclass.DATABLE;
        stepInputData = sid;
    }
    
    @Override
    public void parseStepParameters() throws Exception{
      
    }
    
    
    
    @Override
    public String shortStepDescription(){
      return "Add Flanking Sequence to specified positions";
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "Add Flanking Sequence to specified positions.\n\n"
              + "  1. Input a list of sequences in FASTA format.\n" 
              + "  2. Read the location information that is specified in the header line\n"
              + "  3. add the flanking sequence on either side\n" 
              + "  4. output the new sequence with updated position information\n" 
              + " * \n" 
              + "Header information is in the format:\n" 
              + "  >msy-6088:5|chr4:273774-273790 6088:5|chr4:273774 msy 6088:5|chr4:273774\n" 
              + "  GGTGGAATAGTACCTGT\n" 
              + "i.e. Name|Chromosome:Start or Name|Chromosome:Start-Stop\n" +
" * \n";
    }
    
    
    
    /**
     * This parses out the hashmap containing the run parameters for this step
     * 
     * @param configData
     * @throws Exception 
     */
    @Override
    public void parseConfigurationData(HashMap configData) throws Exception{

        logger.info(STEP_ID_STRING + ": GrabFlankingSequence");
    
        if(configData.get(ID_REF_GENOME)==null) {
            logger.error("<" + ID_REF_GENOME + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_REF_GENOME + "> : Missing Definition in Configuration File");
        }
        this.setReferenceGenome((String) configData.get(ID_REF_GENOME));
        if(this.getReferenceGenome().length() !=3 ){
            logger.error(ID_REF_GENOME + " <" + configData.get(ID_REF_GENOME) + "> must be a 3 letter string");            
            throw new IllegalArgumentException(ID_REF_GENOME + " <" + configData.get(ID_REF_GENOME) + "> must be a 3 letter string");            
        }

        String chk;

        chk = checkParameter("Integer", FIVE_PRIME_FLANK_LEN, Integer.toString((Integer)configData.get(FIVE_PRIME_FLANK_LEN)), "1", "NA", logger);
        if(chk!=null)
            this.setFiveFlankLen((Integer)configData.get(FIVE_PRIME_FLANK_LEN));
        
        chk = checkParameter("Integer", THREE_PRIME_FLANK_LEN, Integer.toString((Integer)configData.get(THREE_PRIME_FLANK_LEN)), "1", "NA", logger);
        if(chk!=null)
            this.setThreeFlankLen((Integer)configData.get(THREE_PRIME_FLANK_LEN));
        
        logger.info("passed");
    }
    
    
    
    
    @Override
    public void execute() throws IOException{
        this.setPaths();
        
        String hostCode = this.getReferenceGenome();
        genomeFasta = new GenomeSeq(hostCode);
        String pathToFasta = getStepInputData().getDataLocations().getGenomeRootFolder()
                + FILESEPARATOR + hostCode + FILESEPARATOR + ReferenceDataLocations.ID_REL_WHOLE_GENSEQ_PATH;
        String genomeFastaFile = this.cleanPath(pathToFasta + FILESEPARATOR + "genome.fa");
        try{
            logger.info("reading genome file <" + genomeFastaFile + ">");
            this.genomeFasta.readFastaGenome(genomeFastaFile);
            logger.info("finished ");
            logger.info("read " + genomeFasta.getNoOfBases() + " bases");
            logger.info("spanning " + genomeFasta.getNoOfChr() + " chromosomes");
        }
        catch(IOException exIO){
            logger.error("exception reading Genome reference file + <" + genomeFastaFile + ">");
            logger.error(exIO.toString());
            throw new IOException(STEP_ID_STRING + ": exception reading Genome reference file + <" + genomeFastaFile + ">");
        }
        
        String gffFeaturesInputFile = "";
        String flankedFastaOutputFile = ""; 
        Boolean f = new File(outFolder).mkdir(); 
        if (f) logger.info("created output folder <" + outFolder + "> for results" );
        
        for(SampleDataEntry sampleData: this.getStepInputData().getSampleData()){
            
            String fastaLine = "";
            try{
                logger.info("reading input features file");
                featureSet = new GFFSet();
                gffFeaturesInputFile = inFolder + FILESEPARATOR + sampleData.getFastqFile1();
                featureSet.readGFF(gffFeaturesInputFile);
                logger.info("-- read " + featureSet.getNoOfEntries() + " entries");
                
                flankedFastaOutputFile =  outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(INFILE_EXTENSION, OUTFILE_EXTENSION);
                BufferedWriter bwFAOut = new BufferedWriter(new FileWriter(new File(flankedFastaOutputFile)));
                    for(GFFEntry feature: featureSet.getGFFEntries()){
                        logger.info(feature.getFeatureID());
                        bwFAOut.write(">" + feature.getFeatureID() 
                                + "|" + feature.getSeqID() + ":" + (feature.getStart()-this.getFiveFlankLen()) 
                                + ":" + feature.getStop()+this.getThreeFlankLen() 
                                + ":" + feature.getStrand() + "\n");
                        bwFAOut.write(this.genomeFasta.getSubSeq(feature.getSeqID(), 
                                feature.getStrand(),
                                feature.getStart()-this.getFiveFlankLen(), 
                                feature.getStop()+this.getThreeFlankLen()) + "\n");
                    }
                bwFAOut.close();
            
            }
            catch(IOException exIO){
                logger.info("error processing Fasta file <" + gffFeaturesInputFile + ">");
                logger.error("error processing Fasta file <" + gffFeaturesInputFile + ">");
                throw new IOException("error processing Fasta file <" + gffFeaturesInputFile + ">\n" + exIO);
             }
        }
        
        logger.info(STEP_ID_STRING + ": completed");
        
    }
    
    
    
    /**
     * parse out a header line for chromosome, strand & position information
     * data is
     * @param headerLine 
     */
    private void parseHeaderLineForLocationData(String headerLine){
        
    }
    
    

    
    /**
     * check data makes sense 
     * @throws IOException, NullPointerException
     * 
     */
    @Override
    public void verifyInputData() throws IOException, NullPointerException{
        
        logger.info(STEP_ID_STRING + ": verify input data");        
        this.setPaths();
                
        
        String pathToFasta = getStepInputData().getDataLocations().getGenomeRootFolder()
                + FILESEPARATOR + this.getReferenceGenome() + FILESEPARATOR + ReferenceDataLocations.ID_REL_WHOLE_GENSEQ_PATH;
        String genomeFastaFile = this.cleanPath(pathToFasta + FILESEPARATOR + "genome.fa");
        if (new File(genomeFastaFile).exists()==false){
            logger.error("no fasta file was found for reference genome <" 
                    + this.getReferenceGenome() + "> at location <" 
                    + genomeFastaFile + ">");
            throw new IOException("no fasta file was found for reference genome <" 
                    + this.getReferenceGenome() + "> at location <" 
                    + genomeFastaFile + ">");
        }



        // check the data files
        for(SampleDataEntry sampleData: this.getStepInputData().getSampleData()){
            
            //Fastq 1
            if (sampleData.getFastqFile1()==null) {
                logger.error(STEP_ID_STRING + ": no Fastq1 file specified");
                throw new IOException(STEP_ID_STRING + ": no Fastq1 file specified");
            }
            String fastqFile1 = inFolder + NGSStep.FILESEPARATOR + sampleData.getFastqFile1().replace(".fastq", INFILE_EXTENSION);
            
            if ((new File(fastqFile1)).exists()==false){
                logger.error(STEP_ID_STRING + ": fastq File1 <" 
                  + fastqFile1 + "> does not exist");
                throw new IOException(STEP_ID_STRING + " : fastq File1 <" 
                  + fastqFile1 + "> does not exist");
            }
            if (fastqFile1.toUpperCase().endsWith(INFILE_EXTENSION.toUpperCase())==false)
            {
                logger.error(STEP_ID_STRING + ": incorrect file extension for input file <" 
                  + fastqFile1 + ">.  \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
                throw new IOException(STEP_ID_STRING + ": incorrect file extension for input file <" 
                  + fastqFile1 + ">.  \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
            }
            
            
            //Fastq 2
            if (sampleData.getFastqFile2()==null) continue;
            String fastqFile2 = inFolder + NGSStep.FILESEPARATOR + sampleData.getFastqFile2().replace(".fastq", INFILE_EXTENSION);
            
            if ((new File(fastqFile2)).exists()==false){
                logger.error(STEP_ID_STRING + " : fastq File2 <" 
                  + fastqFile2 + "> does not exist");
                throw new IOException(STEP_ID_STRING + " : fastq File2 <" 
                  + fastqFile2 + "> does not exist");
            }
            if (fastqFile2.toUpperCase().endsWith(INFILE_EXTENSION.toUpperCase())==false)
            {
                logger.error(STEP_ID_STRING + " : incorrect file extension for fastq file 2 <" 
                  + fastqFile2 + ">. \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
                throw new IOException(STEP_ID_STRING + " : incorrect file extension for fastq file 2 <" 
                  + fastqFile2 + ">. \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
            }
                        
            
        }

    }
    
    @Override
    public void verifyOutputData(){
        
    }

    
    @Override
    public NGSStepSubclass getStepSubclass(){
        return STEP_SUBCLASS;
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
        configData.put(FIVE_PRIME_FLANK_LEN, 25);
        configData.put(THREE_PRIME_FLANK_LEN, 50);
        
        return configData;
    }   
    
    
    
    
    

    /**
     * @return the fiveFlank
     */
    public Integer getFiveFlankLen() {
        return fiveFlankLen;
    }

    /**
     * @param fiveFlank the fiveFlank to set
     */
    public void setFiveFlankLen(Integer fiveFlank) {
        this.fiveFlankLen = fiveFlank;
    }

    /**
     * @return the threeFlank
     */
    public Integer getThreeFlankLen() {
        return threeFlankLen;
    }

    /**
     * @param threeFlank the threeFlank to set
     */
    public void setThreeFlankLen(Integer threeFlank) {
        this.threeFlankLen = threeFlank;
    }

    /**
     * @return the rootDataFolder
     */
    public String getRootDataFolder() {
        return rootDataFolder;
    }

    /**
     * @param rootDataFolder the rootDataFolder to set
     */
    public void setRootDataFolder(String rootDataFolder) {
        this.rootDataFolder = rootDataFolder;
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
