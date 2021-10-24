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
import no.uio.medisin.bag.core.sequence.SimpleSeq;
import no.uio.medisin.bag.ngssmallrna.pipeline.ReferenceDataLocations;
import no.uio.medisin.bag.core.reference.GFFEntry;
import no.uio.medisin.bag.core.reference.GFFSet;
import no.uio.medisin.bag.core.sequence.GenomeSeq;
import no.uio.medisin.bag.core.mapping.MappedRead;
import no.uio.medisin.bag.core.mapping.SAMEntry;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import no.uio.medisin.bag.core.sequence.Strand;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * this step performs additional analyses on a sample set of mapped reads in a
 * SAM file. It can look at features such as
 *
 * 1. distribution of start and end positions 2. distribution of mapped read
 * length
 * 2. overlap of reads with coding regions
 *
 * 
 * However, it only works on the entire dataset, so it can take a very long time
 * to process large read collections
 * 
 * @author sr
 */
public class StepAnalyzeSAMforStartPositions extends NGSStep{
    
    static Logger logger = LogManager.getLogger();
    
    public static final NGSStepSubclass STEP_SUBCLASS           = NGSStepSubclass.DATABLE;

    public  static final String     STEP_ID_STRING      = "AnalyzeSAMforStartPositions";
    private static final String     ID_BLEED            = "bleed";
    private static final String     ID_BASELINE         = "baselinePercent";
    private static final String     ID_MIRBASE_VERSION  = "mirbaseVersion";
    private static final String     ID_REF_GENOME       = "host";
    private static final String     ID_SHORTEST_FEATURE = "shortestFeature";
    private static final String     ID_LONGEST_FEATURE  = "longestFeature";
    private static final String     ID_MIN_COUNTS       = "minCounts";
    private static final String     ID_SEPARATION       = "featureSeparation";
    private static final String     ID_FEATURE_TYPES    = "featureTypes";
    
    
    private static final String     INFILE_EXTENSION    = ".trim.clp.gen.sam";
    private static final String     POS_FILE_EXT        = ".trim.clp.gen.pos.tsv";
    private static final String     FEAT_FILE_EXT       = ".trim.clp.gen.features.tsv";
    
    private static final String     RAW_INPUT_EXTENSION = ".fastq.gz";

    private static final int        COVERAGE_SPAN       = 200;
    
    
    private int                     locationBleed       = 2;
    private int                     baselinePercent     = 5;
    private int                     miRBaseRelease      = 20;
    private String                  ReferenceGenome     = "";
    private int                     shortestRead        = 0;
    private int                     longestRead         = 0;
    private int                     minCounts           = 0;
    private int                     separation          = 0;

    
    private GenomeSeq               genomeFasta;
    private GFFSet                  gffSet              = new GFFSet();    
    private ArrayList<MappedRead>   mappedReads         = new ArrayList<>();
    
    int[]                           coverage5           = new int[COVERAGE_SPAN];    
    int[]                           coverage3           = new int[COVERAGE_SPAN];
    
    private GFFSet                  featureSet          = new GFFSet(); // stores the identified features
    private ArrayList<String>       featureTypes        = new ArrayList<>();
    private ArrayList<String>       featureStrings      = new ArrayList<>();
    

    public StepAnalyzeSAMforStartPositions() {
        classSubtype = NGSStepSubclass.DATABLE;
    }

    /**
     * 
     * @param sid 
     */
    public StepAnalyzeSAMforStartPositions(InputDataForStep sid) {
        classSubtype = NGSStepSubclass.DATABLE;
        stepInputData = sid;
    }

    
    @Override
    public void parseStepParameters() throws Exception{
      
    }
    
    @Override
    public String shortStepDescription(){
      return "investigate read distribution across defined features";
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "investigate read distribution across defined features.\n\n"
              + "The step performs additional analyses on a sample set of mapped reads in a SAM file. "
              + "It can look at features such as\n\n" 
              + "   1. distribution of start and end positions 2. distribution of mapped read length\n" 
              + "   2. overlap of reads with coding regions\n\n" +
                "However, it only works on the entire dataset, so it can take a very long time to process "
              + "large read collections";
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
        
        if(configData.get(ID_BLEED)==null) {
            logger.error("<" + ID_BLEED + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_BLEED + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_SHORTEST_FEATURE)==null) {
            logger.error("<" + ID_SHORTEST_FEATURE + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_SHORTEST_FEATURE + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_LONGEST_FEATURE)==null) {
            logger.error("<" + ID_LONGEST_FEATURE + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_LONGEST_FEATURE + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_MIN_COUNTS)==null) {
            logger.error("<" + ID_MIN_COUNTS + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_MIN_COUNTS + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_BASELINE)==null) {
            logger.error("<" + ID_BASELINE + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_BASELINE + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_MIRBASE_VERSION)==null) {
            logger.error("<" + ID_MIRBASE_VERSION + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_MIRBASE_VERSION + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_REF_GENOME)==null) {
            logger.error("<" + ID_REF_GENOME + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_REF_GENOME + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_SEPARATION)==null) {
            logger.error("<" + ID_SEPARATION + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_SEPARATION + "> : Missing Definition in Configuration File");
        }
        
        if(configData.get(ID_FEATURE_TYPES)!=null)
            this.setFeatureTypes((ArrayList<String>)configData.get(ID_FEATURE_TYPES));


        String chk;
      
        chk = checkParameter("Integer", ID_BASELINE, Integer.toString((Integer)configData.get(ID_BASELINE)), "0", "100", logger);
        if(chk!=null)
            this.setBaselinePercent((Integer)configData.get(ID_BASELINE));

        
        chk = checkParameter("Integer", ID_SHORTEST_FEATURE, Integer.toString((Integer)configData.get(ID_SHORTEST_FEATURE)), "0", "100", logger);
        if(chk!=null)
            this.setShortestRead((Integer)configData.get(ID_SHORTEST_FEATURE));

        

        chk = checkParameter("Integer", ID_LONGEST_FEATURE, Integer.toString((Integer)configData.get(ID_LONGEST_FEATURE)), "0", "NA", logger);
        if(chk!=null)
            this.setLongestRead((Integer)configData.get(ID_LONGEST_FEATURE));


        chk = checkParameter("Integer", ID_MIN_COUNTS, Integer.toString((Integer)configData.get(ID_MIN_COUNTS)), "0", "NA", logger);
        if(chk!=null)
            this.setMinCounts((Integer)configData.get(ID_MIN_COUNTS));


        chk = checkParameter("Integer", ID_BLEED, Integer.toString((Integer)configData.get(ID_BLEED)), "0", "NA", logger);
        if(chk!=null)
            this.setLocationBleed((Integer)configData.get(ID_BLEED));

        
        

        chk = checkParameter("Integer", ID_MIRBASE_VERSION, Integer.toString((Integer)configData.get(ID_MIRBASE_VERSION)), "6", "NA", logger);
        if(chk!=null)
            this.setMiRBaseRelease((Integer)configData.get(ID_MIRBASE_VERSION));
        

        this.setReferenceGenome((String) configData.get(ID_REF_GENOME));
        if(this.getReferenceGenome().length() !=3 ){
            logger.error(ID_REF_GENOME + " <" + configData.get(ID_REF_GENOME) + "> must be a 3 letter string");            
            throw new IllegalArgumentException(ID_REF_GENOME + " <" + configData.get(ID_REF_GENOME) + "> must be a 3 letter string");            
        }

        
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
        File f = new File(pathToAnnotation + FILESEPARATOR + "genes.gtf");
        if (new File(pathToAnnotation + FILESEPARATOR + "genes.gtf").exists()) {
            annotationFile = pathToAnnotation + FILESEPARATOR + "genes.gtf";
        } else if (new File(pathToAnnotation + FILESEPARATOR + "genes.gff").exists()) {
            annotationFile = pathToAnnotation + FILESEPARATOR + "genes.gff";
        }
        
        try {
            gffSet.readGFF(annotationFile);
        } catch (Exception exIO) {
            logger.error("Exception trying to read Annotation file ");
            logger.error(exIO);
            throw new IOException(STEP_ID_STRING + ": exception reading Genome reference file  <" + genomeFastaFile + ">");
        }

        /**
         * read and parse the SAM files
         */
        
        int shortestFeature = this.getShortestRead();
        int longestFeature  = this.getLongestRead();
        
        Boolean fA = new File(outFolder).mkdir();        
        if (fA) {
            logger.info("created output folder <" + outFolder + "> for results");
        }
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
                int featureCount = 0;
        String samLine = null;
        while (itSD.hasNext()) {
            SampleDataEntry sampleData = (SampleDataEntry) itSD.next();
            
            String positionFile = this.cleanPath(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, POS_FILE_EXT));
            String featureOutFile = this.cleanPath(outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, FEAT_FILE_EXT));

            String samInputFile = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION));
            
            logger.info("sam input file is " + samInputFile);
            logger.info("results will be written to " + positionFile);
            
            samLine = null;
            int mapCount = 0;
            try {
                BufferedReader brSAM = new BufferedReader(new FileReader(new File(samInputFile)));
                while ((samLine = brSAM.readLine()) != null) {
                    
                    SAMEntry samEntry = new SAMEntry(samLine);

                    if (samEntry.getqName() == null) {
                        logger.info("broken    "+samLine);
                    }
                    
//                    logger.info("broken    "+samEntry.getqName().split("-").length);
//                    logger.info("broken    "+samEntry.getqName());

                    logger.debug(samLine);
                    if (samEntry.isHeaderLine() == false && samEntry.isMappedRead()) {
//                        logger.info("broken    "+samEntry.getqName());
                        mappedReads.add(new MappedRead(samEntry.getStartPos(), samEntry.getEndPos(),

                                samEntry.getrName(), samEntry.getStrand().toString(), Integer.parseInt(samEntry.getqName().split("-")[1])));

                    }
                    
                }
                brSAM.close();
                logger.debug("read " + mappedReads.size() + " mapped entries");
                
            } catch (IOException smIO) {
                logger.error("error parsing SAM file " + samInputFile);
                logger.error(samLine);
                logger.error(smIO);
                throw new IOException(STEP_ID_STRING + ": error parsing SAM file <" + samInputFile + ">");
            }

            
            featureCount = processSAMReads(featureOutFile, featureCount);
            logger.info("completed parsing");
        }
        
        
        /*
          compare the predicted features to the GFF annotated to see how much 
          overlap we have
        */
        String featVsGFFFile = outFolder + FILESEPARATOR + getStepInputData().getProjectID() + ".gff_vs_feat.tsv";
        try{
            BufferedWriter bwFQ = new BufferedWriter(new FileWriter(new File(featVsGFFFile)));
                for(GFFEntry gffEntry: gffSet.getGFFEntries()){
                    for(GFFEntry feature: featureSet.getGFFEntries()){
                        if(gffEntry.doesFeatureContainRegion(feature)){
                            bwFQ.write("PASS\t" + gffEntry.getType() + "\t" + feature.toGFF3String() + "\n");
                        }else{
                            //bwFQ.write("FAIL\t" + gffEntry.getType() + "\t" + feature.toGFF3String() + "\n");
                        }                            
                    }
                }
            bwFQ.close();
        }
        catch(IOException exIO){
            logger.error("error writing frequency file <" + featVsGFFFile + ">");
            logger.error(exIO);
            throw new IOException(STEP_ID_STRING + ": error writing frequency file <" + featVsGFFFile + ">");
        }
        
        
        String featureFile = outFolder + FILESEPARATOR + getStepInputData().getProjectID() + ".features.tsv";
        try{
            BufferedWriter bwFT = new BufferedWriter(new FileWriter(new File(featureFile)));
                featureSet.writeFeaturesAsGFF3(bwFT, getStepInputData().getProjectID());
            bwFT.close();
        }
        catch(IOException exIO){
            logger.error("error writing feature file <" + featureFile + ">");
            logger.error(exIO);
            throw new IOException(STEP_ID_STRING + ": error writing feature tsv file <" + featureFile + ">");
        }

        String freqFile = outFolder + FILESEPARATOR + getStepInputData().getProjectID() + ".freq.tsv";
        try{
            BufferedWriter bwFQ = new BufferedWriter(new FileWriter(new File(freqFile)));
                featureSet.writeLengthDistribution(bwFQ, 1, longestFeature);
            bwFQ.close();
        }
        catch(IOException exIO){
            logger.error("error writing frequency file <" + freqFile + ">");
            logger.error(exIO);
            throw new IOException(STEP_ID_STRING + ": error writing frequency file <" + freqFile + ">");
        }

        
        
        String fastaFile = outFolder + FILESEPARATOR + getStepInputData().getProjectID() + ".fasta";
        try{
            BufferedWriter bwFA = new BufferedWriter(new FileWriter(new File(fastaFile)));
                featureSet.writeFeaturesAsFastA(bwFA);
            bwFA.close();
        }
        catch(IOException exIO){
            logger.error("error writing fasta feature file <" + fastaFile + ">");
            logger.error(exIO);
            throw new IOException(STEP_ID_STRING + ": error writing fasta feature file <" + fastaFile + ">");
        }
        
        String posFile = outFolder + FILESEPARATOR + getStepInputData().getProjectID() + ".pos.tsv";
        try{
            String headerString = "";
            BufferedWriter bwPS = new BufferedWriter(new FileWriter(new File(posFile)));
                for(String featureType: this.getFeatureTypes()){
                    headerString = headerString.concat(featureType + "\tcounts\t");
                }    
                bwPS.write(headerString.substring(0, headerString.length()-1) + "\n");
                
                for(String featureString: featureStrings){
                    bwPS.write(featureString + "\n");
                }
            bwPS.close();
        }
        catch(IOException exIO){
            logger.error("error writing start position file <" + fastaFile + ">");
            logger.error(exIO);
            throw new IOException(STEP_ID_STRING + ": error writing fasta position file <" + fastaFile + ">");
        }
        
        
        logger.info(STEP_ID_STRING + ": completed");
        
        

    }

    
    
    /**
     *
     * condense overlapping reads into single features, i.e., those
     * which have start/stop positions that vary by <= ´bleed´. Start a
     * new feature when the separation between current read and feature
     * is great than ´separation´.
     *
     * if a read is going to be used to identify a feature, then it is
     * probably sufficient for it to arise in only one sample, because
     * it won´t be identified in a differential expression analysis
     * anyway
     *
     * Also need to accommodate the possibility of two features in the same
     * position but on different strands, therefore we sort by strand and
     * position
     * 
     * @param featureOutFile
     * @param featureCount
     * @return
     * @throws IOException 
     */
    private int processSAMReads(String featureOutFile, int featureCount) throws IOException{
        
        logger.info("processing SAM reads");
        String hostCode = this.getReferenceGenome();
        int shortestFeature = this.getShortestRead();            
        int longestFeature  = this.getLongestRead();
        
        try {
            logger.info("sorting entries...");
            Collections.sort(mappedReads);
            logger.info("...done");

            int bleed = this.getLocationBleed();

            logger.info("preparing feature strings...");
            MappedRead mappedRead = ((MappedRead) mappedReads.get(0));
            //logger.debug(mappedRead.toString());
            String currentChr = mappedRead.getChr();
            Strand currentStrand = mappedRead.getStrand();
            String featureString = "";
            for(String featureType: this.getFeatureTypes()){
                featureString = featureString.concat(findOverlappingFeature(mappedRead, featureType) + "\t");                    
            }
            featureStrings.add(featureString.substring(0, featureString.length()-1));
            logger.info("...done");

            /*
             need to accommodate the possibility of two features in the
             same position but on different strands
             */
            int currentStart5 = 100000000;
            int currentStop5 = -1;
            int currentStart3 = -1;
            int currentStop3 = 100000000;

            /*
             using an int array might not be the best way to proceed, but try for now
             */
            int coverage5Start = -1;
            int coverage3Start = -1;

            if (currentStrand.equals(Strand.PLUS)) {
                currentStart5 = mappedRead.getStartPos();
                currentStop5 = mappedRead.getEndPos();
                coverage5Start = currentStart5 - COVERAGE_SPAN / 2;
            } else {
                currentStart3 = mappedRead.getStartPos();
                currentStop3 = mappedRead.getEndPos();                    
                coverage3Start = currentStart3 - COVERAGE_SPAN / 2;
            }


            BufferedWriter bwFT = new BufferedWriter(new FileWriter(new File(featureOutFile)));
            bwFT.write("ID\tChr\t \tStart\tStop\tLength\tCoverage\tDispersion\n");

            logger.info("processing mapped reads...");
            Iterator itMR = mappedReads.iterator();
            Boolean startNewFeature3 = false;
            Boolean startNewFeature5 = false;
            while (itMR.hasNext()) {
                mappedRead = (MappedRead) itMR.next();
                logger.debug(mappedRead);
                featureString = "";
                for(String featureType: this.getFeatureTypes()){
                    featureString = featureString.concat(findOverlappingFeature(mappedRead, featureType) + "\t");                    
                }
                featureStrings.add(featureString.substring(0, featureString.length()-1));
                //logger.debug(mappedRead.toString());
                if (mappedRead.getChr().equals(currentChr) == false || mappedRead.getStrand().equals(currentStrand) == false) {
                    if (mappedRead.getChr().equals(currentChr) == false) {
                        //logger.debug("different chr. start new feature");
                    }
                    if (mappedRead.getStrand().equals(currentStrand) == false) {
                        logger.debug("different strand. start new feature");                            
                    }
                    currentStrand = mappedRead.getStrand();
                    currentChr = mappedRead.getChr();
                    if (currentStrand == Strand.PLUS) {
                        startNewFeature3 = false;
                        startNewFeature5 = false;
                        currentStart5 = mappedRead.getStartPos();
                        currentStop5 = mappedRead.getEndPos();
                        coverage5Start = currentStart5 - COVERAGE_SPAN / 2;
                        addCounts5(currentStart5 - coverage5Start, mappedRead.getEndPos() - coverage5Start, mappedRead.getCount());
                    } else {
                        startNewFeature3 = false;
                        startNewFeature5 = false;
                        currentStart3 = mappedRead.getStartPos();
                        currentStop3 = mappedRead.getEndPos();
                        coverage3Start = currentStart3 + COVERAGE_SPAN / 2;
                        addCounts3(coverage3Start - mappedRead.getStartPos(), coverage3Start - mappedRead.getEndPos(), mappedRead.getCount());
                    }

                } else {
                    switch (currentStrand) {
                        case PLUS:
                            if (mappedRead.getStartPos() - currentStart5 > separation) {
                                //logger.debug("read too far away. start new 5' feature");                                    
                                startNewFeature5 = true;
                            } else {
                                if (mappedRead.getStartPos() < currentStart5) {
                                    //logger.debug("extend 5´ start from " + currentStart5 + " to " + mappedRead.getStartPos());
                                    currentStart5 = mappedRead.getStartPos();
                                    if( currentStart5 <= coverage5Start)
                                        startNewFeature5 = true;
                                }
                                if (mappedRead.getEndPos() > currentStop5) {
                                    //logger.debug("extend 5´ end from " + currentStop5 + " to " + mappedRead.getEndPos());
                                    currentStop5 = mappedRead.getEndPos();
                                    if( currentStop5 >= coverage5Start + COVERAGE_SPAN )
                                        startNewFeature5 = true;
                                }
                                try {
                                    //logger.debug(printCoverage5(currentStart5 - coverage5Start, mappedRead.getEndPos() - coverage5Start, "|"));
                                    this.addCounts5(currentStart5 - coverage5Start, mappedRead.getEndPos() - coverage5Start, mappedRead.getCount());
                                    //logger.debug(printCoverage5(currentStart5 - coverage5Start, mappedRead.getEndPos() - coverage5Start, "|"));
                                } catch (ArrayIndexOutOfBoundsException exAB) {
                                    logger.error("5' Array out of bounds");
                                    throw new IOException("5' Array out of bounds");
                                }
                            }
                            break;


                        case MINUS:
                            if ((mappedRead.getStartPos() - currentStart3 > separation) || (currentStop3 - mappedRead.getEndPos() > separation)) {
                                //logger.debug("read too far away. start new 3' feature");
                                startNewFeature3 = true;
                            } else {
                                //logger.debug("add read to existing feature");
                                if (mappedRead.getStartPos() > currentStart3) {
                                    //logger.debug("extend 3´ start from " + currentStart3 + " to " + mappedRead.getStartPos());
                                    currentStart3 = mappedRead.getStartPos();
                                    if(currentStart3 <= coverage3Start)
                                        startNewFeature3 = true;
                                }
                                if (mappedRead.getEndPos() < currentStop3) {
                                    //logger.debug("extend 3´ end from " + currentStop3 + " to " + mappedRead.getEndPos());
                                    currentStop3 = mappedRead.getEndPos();
                                    if(currentStop3 >= coverage3Start + COVERAGE_SPAN) 
                                        startNewFeature3 = true;
                                }
                                try {
                                    //logger.debug(printCoverage3(coverage3Start - currentStop3, coverage3Start - currentStart3, "|"));
                                    this.addCounts3(coverage3Start - mappedRead.getStartPos(), coverage3Start - mappedRead.getEndPos(), mappedRead.getCount());
                                    //logger.debug(printCoverage3(coverage3Start - currentStop3, coverage3Start - currentStart3, "|"));
                                } catch (ArrayIndexOutOfBoundsException exAB) {
                                    logger.error("3' Array out of bounds");
                                    logger.error(mappedRead.toString());
                                    logger.error(exAB);
                                    throw new IOException("3' Array out of bounds" + "\n" + mappedRead.toString());
                                }
                            }
                            break;


                        case UNKNOWN:
                            logger.error("no Strand information for read, cannot process ");
                            logger.error(mappedRead.toString());
                            break;

                    }

                }

                if (startNewFeature5) {
                    if(currentStop5 - currentStart5 + 1 <= longestFeature 
                            && currentStop5 - currentStart5 + 1>= shortestFeature 
                            && countCoverage5(currentStart5 - coverage5Start, currentStop5 - coverage5Start) >= minCounts){
                        bwFT.write(featureCount + "\t" + currentChr + "\t" + currentStrand + "\t" + currentStart5 + "\t" + currentStop5 + "\t"
                                + (currentStop5 - currentStart5 + 1) + "\t" + this.countCoverage5(currentStart5 - coverage5Start, currentStop5 - coverage5Start)
                                + "\t" + this.countDispersion5(currentStart5 - coverage5Start, currentStop5 - coverage5Start) + "\n");
                        if(featureSet.doesRegionContainFeature(currentStart5, currentStop5, currentStrand, currentChr, bleed)==false){
                            GFFEntry newEntry = new GFFEntry(
                                    currentChr,
                                    "srp",
                                    "smallRNA",
                                    currentStart5, 
                                    currentStop5,
                                    ".",
                                    currentStrand.toString(),
                                    "0",
                                    "ID=" + hostCode  + "-" + Integer.toString(featureCount) + ":5|" + currentChr  + ":" + currentStart5+ "-" + currentStop5
                                    );
                            String featureSeq = genomeFasta.getSubSeq(currentChr, currentStrand, currentStart5, currentStop5);
                            newEntry.setSequence(featureSeq);
                            featureSet.addEntry(newEntry);
                        }                            
                    }

                    currentStart5 = mappedRead.getStartPos();
                    currentStop5 = mappedRead.getEndPos();

                    coverage5Start = currentStart5 - COVERAGE_SPAN / 2;
                    clearCoverage5();
                    try {
                        this.addCounts5(mappedRead.getStartPos() - coverage5Start, mappedRead.getEndPos() - coverage5Start, mappedRead.getCount());
                    } catch (ArrayIndexOutOfBoundsException exAB) {
                        logger.error("5' Array out of bounds when writing data");
                        throw new IOException(STEP_ID_STRING + ": 5' Array out of bounds when writing data");
                    }

                    currentStrand = mappedRead.getStrand();
                    currentChr = mappedRead.getChr();
                    featureCount++;
                    startNewFeature5 = false;

                } else {                                         
                    if (startNewFeature3 ) {
                        if(currentStop3 - currentStart3 + 1<= longestFeature 
                            && currentStop3 - currentStart3 + 1>= shortestFeature
                            && countCoverage3(coverage3Start - currentStop3, coverage3Start - currentStart3) > minCounts){
                            bwFT.write(featureCount + "\t" + currentChr + "\t" + currentStrand + "\t" + currentStart3 + "\t" + currentStop3 + "\t"
                                    + (currentStop3 - currentStart3 + 1) + "\t" + this.countCoverage3(coverage3Start - currentStop3, coverage3Start - currentStart3)
                                    + "\t" + this.countDispersion3(coverage3Start - currentStop3, coverage3Start - currentStart3) + "\n");
                            if(featureSet.doesRegionContainFeature(currentStart3, currentStop3, currentStrand, currentChr, bleed)==false){
     //                           GFFEntry newEntry = new GFFEntry(hostCode + "-" + Integer.toString(featureCount) + ":3|" + currentChr  + ":" + currentStart5+ "-" + currentStop5, currentStrand.toString(), currentChr, currentStart3, currentStop3);
                                GFFEntry newEntry = new GFFEntry(
                                        currentChr,
                                        "srp",
                                        "smallRNA",
                                        currentStart5, 
                                        currentStop5,
                                        ".",
                                        currentStrand.toString(),
                                        "0",
                                        "ID=" + hostCode + "-" + Integer.toString(featureCount) + ":3|" + currentChr  + ":" + currentStart5+ "-" + currentStop5
                                        );
                                String featureSeq = SimpleSeq.complement(genomeFasta.getSubSeq(currentChr, currentStrand, currentStart3, currentStop3));
                                newEntry.setSequence(featureSeq);
                                featureSet.addEntry(newEntry);                                    
                            }

                        }

                        currentStart3 = mappedRead.getStartPos();
                        currentStop3 = mappedRead.getEndPos();

                        coverage3Start = currentStart3 + COVERAGE_SPAN / 2;
                        clearCoverage3();
                        try {
                            this.addCounts3(coverage3Start - mappedRead.getStartPos(), coverage3Start - mappedRead.getEndPos(), mappedRead.getCount());
                        } catch (ArrayIndexOutOfBoundsException exAB) {
                            logger.error("3' Array out of bounds when writing data");
                            throw new IOException(STEP_ID_STRING + ": 3' Array out of bounds when writing data");
                        }
                        currentStrand = mappedRead.getStrand();
                        currentChr = mappedRead.getChr();                           
                        featureCount++;
                        startNewFeature3 = false;
                    }
                }

            }
            bwFT.close();
            logger.info("...done");
        } catch (IOException exIO) {
            logger.error(exIO);
            throw new IOException(STEP_ID_STRING + ": error writing feature details file <" + featureOutFile + ">");
        }
        return featureCount;
    }
    
    
    
    
    /**
     * Is there a feature in the GFFEntry table that overlaps this read?
     * 
     * @param queryRead
     * @param featureType
     * @return 
     */
    private int checkFeatureOverlap(MappedRead queryRead, String featureType, int bleed){
        GFFEntry gffEntry = gffSet.findOverlappingFeature(queryRead, featureType, bleed);
        if (gffEntry!=null)
            return gffEntry.getStart() - queryRead.getStartPos();
        else
            return -9999;
    }
    
    
    



    /**
     * Is there a feature in the GFFEntry table that overlaps this read?
     * 
     * @param queryRead
     * @param featureType
     * @return 
     */
    private String findOverlappingFeature(MappedRead queryRead, String featureType){
        GFFEntry gffEntry = gffSet.doesFeatureContainRegion(queryRead, featureType);
        if (gffEntry!=null)
            return Integer.toString(gffEntry.getStart() - queryRead.getStartPos()) + "\t" + queryRead.getCount();
        else
            return "-8888\t" + queryRead.getCount();
    }
    
    
    
    
    
    
    /**
     * calculate average count across the mapped region for + Strand
     *
     * @param start
     * @param stop
     * @return
     */
    private double countCoverage5(int start, int stop) {
        int countTotal = 0;
        for (int i = start; i < stop; i++) {
            countTotal += coverage5[i];
        }
        return (double) countTotal / (double) (stop - start + 1);
    }

    /**
     * add counts to 3' coverage
     *
     * @param start
     * @param stop
     * @param counts
     */
    private void addCounts3(int start, int stop, int counts) {
        
        for (int i = stop; i < start; i++) {
            coverage3[i] += counts;
        }
        
    }

    /**
     * add counts to 5' coverage
     *
     * @param start
     * @param stop
     * @param counts
     */
    private void addCounts5(int start, int stop, int counts) {
        
        for (int i = start; i < stop; i++) {
            coverage5[i] += counts;
        }
        
    }

    /**
     * calculate average count across the mapped region for - Strand
     *
     * @param start
     * @param stop
     * @return
     */
    private double countCoverage3(int start, int stop) {
        int countTotal = 0;
        for (int i = start; i < stop; i++) {
            countTotal += coverage3[i];
        }
        return (double) countTotal / (double) (stop - start + 1);
    }

    /**
     * return string showing coverage over specified 3' region
     *
     * @param start
     * @param stop
     * @return
     */
    private String printCoverage3(int start, int stop, String delimiter) {
        String coverageStr = "";
        for (int i = start; i < stop; i++) {
            coverageStr = coverageStr.concat((Integer.toString(coverage3[i])) + delimiter);
        }
        return coverageStr;
    }

    /**
     * return string showing coverage over specified 3' region
     *
     * @param start
     * @param stop
     * @return
     */
    private String printCoverage5(int start, int stop, String delimiter) {
        String coverageStr = "";
        for (int i = start; i < stop; i++) {
            coverageStr = coverageStr.concat((Integer.toString(coverage5[i])) + delimiter);
        }
        return coverageStr;
    }

    /**
     * reset the 3' coverage array
     */
    private void clearCoverage3() {
        for (int i = 0; i < coverage3.length; i++) {
            coverage3[i] = 0;
        }        
    }

    /**
     * estimate dispersion of counts across the feature for + strand
     *
     * @param start
     * @param stop
     * @return
     */
    private double countDispersion5(int start, int stop) {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (int i = start; i < stop; i++) {
            stats.addValue(coverage5[i]);
        }
        return stats.getStandardDeviation();
    }

    /**
     * reset the 5' coverage array
     */
    private void clearCoverage5() {
        for (int i = 0; i < coverage5.length; i++) {
            coverage5[i] = 0;
        }        
    }

    /**
     * estimate dispersion of counts across the feature for - strand
     *
     * @param start
     * @param stop
     * @return
     */
    private double countDispersion3(int start, int stop) {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (int i = start; i < stop; i++) {
            stats.addValue(coverage3[i]);
        }
        return stats.getStandardDeviation();
    }

    /**
     * Verify Input Data for parsing SAM file for miRNAs
     * @throws IOException
     */    
    @Override
    public void verifyInputData()  throws IOException{

        logger.info(STEP_ID_STRING + " :verify input data");        
        this.setPaths();
        
        
        String gffFileMirBase = this.cleanPath(getStepInputData().getDataLocations().getMirbaseFolder() 
                + FILESEPARATOR + this.getMiRBaseRelease() + FILESEPARATOR + this.getReferenceGenome() + ".gff3");
        if (new File(gffFileMirBase).exists()==false){
            logger.error("no annotation file was found for mirBase HOST:<" 
                    + this.getReferenceGenome() + "> VERSION: <"+ this.getMiRBaseRelease() + "> at location <" 
                    + gffFileMirBase + ">");
            throw new IOException("no annotation file was found for mirBase HOST:<" 
                    + this.getReferenceGenome() + "> VERSION: <"+ this.getMiRBaseRelease() + "> at location <" 
                    + gffFileMirBase + ">");
        }
                
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
        configData.put(ID_MIRBASE_VERSION, 20);
        configData.put(ID_BLEED, 2);
        configData.put(ID_BASELINE, 5);
        configData.put(ID_SHORTEST_FEATURE, 2);
        configData.put(ID_LONGEST_FEATURE, 200);
        configData.put(ID_MIN_COUNTS, 1000);
        configData.put(ID_SEPARATION, 10);
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
     * @return the locationBleed
     */
    public int getLocationBleed() {
        return locationBleed;
    }

    /**
     * @param locationBleed the locationBleed to set
     */
    public void setLocationBleed(int locationBleed) {
        this.locationBleed = locationBleed;
    }

    /**
     * @return the baselinePercent
     */
    public int getBaselinePercent() {
        return baselinePercent;
    }

    /**
     * @param baselinePercent the baselinePercent to set
     */
    public void setBaselinePercent(int baselinePercent) {
        this.baselinePercent = baselinePercent;
    }

    /**
     * @return the miRBaseRelease
     */
    public int getMiRBaseRelease() {
        return miRBaseRelease;
    }

    /**
     * @param miRBaseRelease the miRBaseRelease to set
     */
    public void setMiRBaseRelease(int miRBaseRelease) {
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
