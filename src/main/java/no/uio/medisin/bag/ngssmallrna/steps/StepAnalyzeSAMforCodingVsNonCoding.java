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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import no.uio.medisin.bag.ngssmallrna.pipeline.ReferenceDataLocations;
import no.uio.medisin.bag.core.reference.GFFEntry;
import no.uio.medisin.bag.core.reference.GFFSet;
import no.uio.medisin.bag.core.sequence.GenomeSeq;
import no.uio.medisin.bag.core.mapping.MappedRead;
import no.uio.medisin.bag.core.mapping.SAMEntry;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * this step performs additional analyses on a sample set of mapped reads in a
 * SAM file. 
 *
 * here, it specifically finds the percentage of reads that map to coding
 * and non-coding regions of the genome
 * 
 * This could be done within the @see no.uio.medisin.bag.ngssmallrna.steps.StepAnalyzeSAMforStartPositions
 * class, but there is significant overhead associated with the analysis of the
 * start positions, so we split the non-coding/coding analysis into this separate
 * class
 * 
 * 
 * @author sr
 */
public class StepAnalyzeSAMforCodingVsNonCoding extends NGSStep{
    
    static Logger logger = LogManager.getLogger();
    
    public static final NGSStepSubclass STEP_SUBCLASS           = NGSStepSubclass.DATABLE;

    public  static final String     STEP_ID_STRING      = "StepAnalyzeSAMforCodingVsNonCoding";
    private static final String     ID_REF_GENOME       = "host";
    private static final String     ID_FEATURE_TYPES    = "featureTypes";
    
    
    private static final String     INFILE_EXTENSION    = ".trim.clp.gen.sam";
    private static final String     COUNT_SUMMARY_EXT   = ".trim.clp.gen.count_summary.tsv";
    
    private static final String     RAW_INPUT_EXTENSION = ".fastq.gz";

    
    
    private String                  ReferenceGenome     = "";
    private int                     shortestRead        = 0;
    private int                     longestRead         = 0;
    private int                     minCounts           = 0;
    private int                     separation          = 0;
    private long                    codingReadCount     = 0;
    private long                    nonCodingReadCount  = 0;
    private long                    totalReadCount      = 0;

    
    private GenomeSeq               genomeFasta;
    private GFFSet                  gffSet              = new GFFSet();    
    private ArrayList<MappedRead>   mappedReads         = new ArrayList<>();
    
    
    private ArrayList<String>       featureTypes        = new ArrayList<>();
    

    @Override
    public void parseStepParameters() throws Exception{
      
    }
    
    public StepAnalyzeSAMforCodingVsNonCoding() {
        classSubtype = NGSStepSubclass.DATABLE;
    }

    /**
     * 
     * @param sid 
     */
    public StepAnalyzeSAMforCodingVsNonCoding(InputDataForStep sid) {
        classSubtype = NGSStepSubclass.DATABLE;
        stepInputData = sid;
    }

    
    

    @Override
    public String shortStepDescription(){
      return "finds the percentage of reads that map to coding and non-coding regions of the genome";
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "finds the percentage of reads that map to coding and non-coding regions of the genome.\n\n"
              + "here, it simply finds the percentage of reads that map to coding and non-coding regions "
              + "of the genome. This could be done within the StepAnalyzeSAMforStartPositions class but "
              + "there is significant overhead associated with the analysis of the start positions, "
              + "so we split the non-coding/coding analysis into this separate class";
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


        this.setReferenceGenome((String) configData.get(ID_REF_GENOME));
        if(this.getReferenceGenome().length() !=3 ){
            logger.error(ID_REF_GENOME + " <" + configData.get(ID_REF_GENOME) + "> must be a 3 letter string");            
            throw new IllegalArgumentException(ID_REF_GENOME + " <" + configData.get(ID_REF_GENOME) + "> must be a 3 letter string");            
        }

        if(configData.get(ID_FEATURE_TYPES)!=null)
            this.setFeatureTypes((ArrayList<String>)configData.get(ID_FEATURE_TYPES));
        
        logger.info("passed");
    }
    
        
    
    
    
    /**
     * parse out SAM file to retrieve start and stop information for each read
     *
     *
     * @param filename
     * @throws IOException
     */
    @Override
    public void execute() throws IOException{
        
        logger.info(STEP_ID_STRING + ": execute");                
        this.setPaths();
    
        
        /**
         * read genome fasta
         */
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
        
        /**
         * read genes.gtf or genes.gff3 file
         */
        String annotationFile = "";
        String pathToAnnotation = getStepInputData().getDataLocations().getGenomeRootFolder()
                + FILESEPARATOR + hostCode + ReferenceDataLocations.ID_GENE_ANNOTATION;

        if (new File(pathToAnnotation + FILESEPARATOR + "genes.gtf").exists()) {
            annotationFile = pathToAnnotation + FILESEPARATOR + "genes.gtf";
        } else if (new File(pathToAnnotation + FILESEPARATOR + "genes.gff").exists()) {
            annotationFile = pathToAnnotation + FILESEPARATOR + "genes.gff";
        }
        
        try {
            gffSet.readGFF(annotationFile, this.getFeatureTypes());
        } catch (Exception exIO) {
            logger.error("Exception trying to read Annotation file ");
            logger.error(exIO);
            throw new IOException(STEP_ID_STRING + ": exception reading Genome reference file  <" + genomeFastaFile + ">");
        }

        Boolean fA = new File(outFolder).mkdir();        
        if (fA) {
            logger.info("created output folder <" + outFolder + "> for results");
        }
        
        /**
         * read and parse the SAM files
         */
        String featureOutFile = this.cleanPath(outFolder + FILESEPARATOR + this.getStepInputData().getProjectID() + ".mapping_summary.tsv");
        logger.info("results will be written to " + featureOutFile);
        try {
            BufferedWriter bwFT = new BufferedWriter(new FileWriter(new File(featureOutFile)));
            bwFT.write("ID\tTotal Counts\tMapped Counts\tCoding Counts\tNonCoding Counts\n");

            Iterator itSD = this.getStepInputData().getSampleData().iterator();
            String samLine = null;
            while (itSD.hasNext()) {
                SampleDataEntry sampleData = (SampleDataEntry) itSD.next();


                String samInputFile = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION));

                logger.info("sam input file is " + samInputFile);

                samLine = null;
                totalReadCount = 0;
                codingReadCount = 0;
                nonCodingReadCount = 0;
                try {
                    BufferedReader brSAM = new BufferedReader(new FileReader(new File(samInputFile)));
                    while ((samLine = brSAM.readLine()) != null) {
                        logger.debug(samLine);

                        SAMEntry samEntry = new SAMEntry(samLine);
                        if (samEntry.isHeaderLine()) 
                            continue;
                        totalReadCount += Integer.parseInt(samEntry.getqName().split("-")[1]);
                        if (samEntry.isHeaderLine() == false && samEntry.isMappedRead()) {
                            mappedReads.add(new MappedRead(samEntry.getStartPos(), samEntry.getEndPos(),
                                    samEntry.getrName(), samEntry.getStrand().toString(), Integer.parseInt(samEntry.getqName().split("-")[1])));
                        }

                    }
                    brSAM.close();
                    logger.debug("read " + mappedReads.size() + " mapped entries");
                    logger.info("total mapped reads:\t" + totalReadCount);

                } catch (IOException smIO) {
                    logger.error("error parsing SAM file " + samInputFile);
                    logger.error(samLine);
                    logger.error(smIO);
                    throw new IOException(STEP_ID_STRING + ": error parsing SAM file <" + samInputFile + ">");
                }


                logger.info("sorting entries...");
                Collections.sort(mappedReads);
                Collections.sort(this.gffSet.getGFFEntries());
                logger.info("...done");

                logger.info("processing mapped reads...");

                MappedRead mappedRead;

                int i=0;
                Iterator itMR = mappedReads.iterator();
                while (itMR.hasNext()) {
                    mappedRead = (MappedRead) itMR.next();
                    logger.debug(mappedRead);
                    for(String featureType: this.getFeatureTypes()){
                        if(this.doesReadOverlapFeature(mappedRead, featureType, 0))
                            codingReadCount += mappedRead.getCount();
                        else
                            nonCodingReadCount += mappedRead.getCount();
                    }
                    i++;
                    if(i%1000==0)
                        logger.info(".");
                }
                bwFT.write(sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, COUNT_SUMMARY_EXT) + "\t" 
                        + totalReadCount + "\t" + (codingReadCount + nonCodingReadCount) + "\t" 
                        + codingReadCount + "\t" + nonCodingReadCount + "\n");
                logger.info("...done");
                logger.info("completed parsing");
            }        
        
            bwFT.close();
        } catch (IOException exIO) {
            logger.error(exIO);
            throw new IOException(STEP_ID_STRING + ": error writing feature details file <" + featureOutFile + ">");
        }        
        logger.info(STEP_ID_STRING + ": completed");
        
    }

    
    
    /**
     * check whether the supplied read overlaps a feature of the specified type
     * 
     * @param queryRead
     * @param featureType
     * @param bleed
     * @return 
     */
    private Boolean doesReadOverlapFeature(MappedRead queryRead, String featureType, int bleed){
        
        GFFEntry gffEntry = gffSet.findOverlappingFeature(queryRead, featureType, bleed);
        
        return gffEntry!=null;
        
    }
    
    
    
    /**
     * Verify Input Data for parsing SAM file for miRNAs
     * @throws IOException
     */    
    @Override
    public void verifyInputData()  throws IOException{

        logger.info(STEP_ID_STRING + " :verify input data");        
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
        
        String annotationFile = "";
        String pathToAnnotation = getStepInputData().getDataLocations().getGenomeRootFolder()
                + FILESEPARATOR + this.getReferenceGenome() + ReferenceDataLocations.ID_GENE_ANNOTATION;
        File f = new File(pathToAnnotation + FILESEPARATOR + "genes.gtf");
        if (new File(pathToAnnotation + FILESEPARATOR + "genes.gtf").exists()) {
            annotationFile = pathToAnnotation + FILESEPARATOR + "genes.gtf";
        } else if (new File(pathToAnnotation + FILESEPARATOR + "genes.gff").exists()) {
            annotationFile = pathToAnnotation + FILESEPARATOR + "genes.gff";
        }
        if (annotationFile == null) {
            logger.error("no annotation file was found for reference genome <" 
                    + pathToAnnotation + "> at location <"
                    + annotationFile + ">");
            throw new IOException("no annotation file was found for reference genome <" 
                    + pathToAnnotation + "> at location <"
                    + annotationFile + ">");
        }
        
        
                    
        // check the data files
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            
            //Fastq 1 SAM file
            if (sampleData.getFastqFile1()==null) throw new IOException("no Fastq1 file specified");
            String samFileIn = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION));
            
            if ((new File(samFileIn)).exists()==false){
                logger.error("SAM input file <" 
                  + samFileIn + "> does not exist");
                throw new IOException("SAM input file  <" 
                  + samFileIn + "> does not exist");
            }
            if (samFileIn.toUpperCase().endsWith(INFILE_EXTENSION.toUpperCase())==false)
            {
                logger.error("incorrect file extension for SAM input file <" 
                  + samFileIn + ">.  \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
                throw new IOException("incorrect file extension for SAM input file <" 
                  + samFileIn + ">.  \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
            }
        }
        logger.info("passed");
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
        configData.put(ID_FEATURE_TYPES, new ArrayList<>(Arrays.asList("mRNA", "CDS", "exon")));
        
        return configData;
    }





    @Override
    public void verifyOutputData() {
        
    }

    @Override
    public NGSStepSubclass getStepSubclass(){
        return STEP_SUBCLASS;
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

    /**
     * @return the shortestRead
     */
    public int getShortestRead() {
        return shortestRead;
    }

    /**
     * @param shortestRead the shortestRead to set
     */
    public void setShortestRead(int shortestRead) {
        this.shortestRead = shortestRead;
    }

    /**
     * @return the longestRead
     */
    public int getLongestRead() {
        return longestRead;
    }

    /**
     * @param longestRead the longestRead to set
     */
    public void setLongestRead(int longestRead) {
        this.longestRead = longestRead;
    }

    /**
     * @return the min_counts
     */
    public int getMinCounts() {
        return minCounts;
    }

    /**
     * @param min_counts the min_counts to set
     */
    public void setMinCounts(int min_counts) {
        this.minCounts = min_counts;
    }

    /**
     * @return the separation
     */
    public int getSeparation() {
        return separation;
    }

    /**
     * @param separation the separation to set
     */
    public void setSeparation(int separation) {
        this.separation = separation;
    }

    /**
     * @return the featureTypes
     */
    public ArrayList<String> getFeatureTypes() {
        return featureTypes;
    }

    /**
     * @param featureTypes the featureTypes to set
     */
    public void setFeatureTypes(ArrayList<String> featureTypes) {
        this.featureTypes = featureTypes;
    }


    
}
