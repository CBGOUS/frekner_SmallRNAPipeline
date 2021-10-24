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
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import no.uio.medisin.bag.core.reference.GFFEntry;
import no.uio.medisin.bag.core.reference.GFFSet;
import no.uio.medisin.bag.core.mirna.MiRNASet;
import no.uio.medisin.bag.ngssmallrna.pipeline.ReferenceDataLocations;
import no.uio.medisin.bag.core.mapping.SAMEntry;
import no.uio.medisin.bag.core.targetscan.TargetScanMirFamilyList;
import no.uio.medisin.bag.core.targetscan.TargetScanPredictedTargetList;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;

import org.apache.logging.log4j.Logger;




/**
 * match two smallRNAs by comparing seed regions 
 * The data consists of Query and Reference sets. 
 * 
 * The Query data comprises
 *   a Features file in GFF3 format with location information (start, stop, chromosome and strand)
 *   a FASTA file with sequence information
 *   a Counts file with count information for each features across a set of samples
 * 
 * The Reference data comprises
 *   a Features file in GFF3 format with location information for smallRNAs (i.e. miRNAs)
 *   a TargetScan 'miR Family' file with miRNA family and seed region information
 *   a TargetScan 'predicted targets' file with target information
 * 
 * The execute step consolidates the Query GFF, Fasta & Count files so that each feature can be filtered
 * according to counts, and a seed region can be extracted from the sequence. This seed region can 
 * then be compared to seed regions in the Reference genome, as identified by TargetScan
 * Target information (in terms of Gene and TranscriptIDs) can then be obtained for overlapping seed regions
 * 
 * Note: for now, the seed regions are only matched bases on nt 2 to 8, nothing more complex is performed in
 * terms of matching.
 * 
 * @author sr
 */

public class StepMatchSmallRNAsBySeedRegions extends NGSStep{
    
    static Logger                       logger = LogManager.getLogger();
    
    public static final NGSStepSubclass STEP_SUBCLASS           = NGSStepSubclass.DATABLE;

    public static final String          STEP_ID_STRING          = "MatchSmallRNAsBySeedRegions";
    private static final String         ID_MIRBASE_VERSION      = "mirbaseVersion";
    private static final String         ID_REF_GENOME           = "host";
    private static final String         ID_QUERY_GENOME         = "query";    
    private static final String         ID_THREADS              = "noOfThreads";
    private static final String         ID_MINCOUNTS            = "minCounts";
    
    private static final String         INFILE_EXTENSION        = ".fastq.gz";
    private static final String         OUTFILE_EXTENSION       = ".fastq";
    
    private static final String         QUERYFA_EXTENSION       = ".fasta";
    private static final String         QUERYCOUNTS_EXTENSION   = ".merged.mirna_counts.tsv";
    private static final String         QUERYFEATURES_EXTENSION = ".features.tsv";
    
    MiRNASet                            mirBaseSet              = new MiRNASet(); 
    TargetScanMirFamilyList             tScanMirFamilies        = new TargetScanMirFamilyList();
    TargetScanPredictedTargetList       tscanPredictedTargets   = new TargetScanPredictedTargetList();
    
    private int                         noOfThreads             = 4;
    private int                         minCounts               = 1000;
    private String                      unzipSoftware           = "";
    private int                         miRBaseRelease          = 20;
    private String                      referenceGenome         = "";
    
    
    private String                      queryGenome             = "";
    private String                      queryCountDataFile      = "";
    private String                      queryFastaFile          = "";
    private String                      queryFeatureFile        = "";
    private String                      hostVersusQuerySAMFile  = "";
    
    private ArrayList<String>           querySeedSet;
    private GFFSet                      queryFeatures;
    
    
    

    
    public StepMatchSmallRNAsBySeedRegions(){
        classSubtype = NGSStepSubclass.DATABLE;
    }
    

    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepMatchSmallRNAsBySeedRegions(InputDataForStep sid){
        classSubtype = NGSStepSubclass.DATABLE;
       stepInputData = sid;
    }
    
    @Override
    public void parseStepParameters() throws Exception{
      
    }
    


    @Override
    public String shortStepDescription(){
      return "match two smallRNAs by comparing seed regions";
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "match two smallRNAs by comparing seed regions \n\n" +
        " The data consist of Query and Reference sets. \n" +
        " \n" +
        " The Query data comprises:\n\n" +
        "  - a Features file in GFF3 format with location information (start, stop, chromosome and strand)\n" +
        "  - a FASTA file with sequence information\n" +
        "  - a Counts file with count information for each features across a set of samples\n" +
        "   \n" +
        "  The Reference data comprises:\n\n" +
        "   - a Features file in GFF3 format with location information for smallRNAs (i.e. miRNAs)\n" +
        "   - a TargetScan 'miR Family' file with miRNA family and seed region information\n" +
        "   - a TargetScan 'predicted targets' file with target information\n" +
        " \n" +
        " The execute step consolidates the Query GFF, Fasta & Count files so that each feature can be filtered\n" +
        " according to counts, and a seed region can be extracted from the sequence. This seed region can \n" +
        " then be compared to seed regions in the Reference genome, as identified by TargetScan\n" +
        " Target information (in terms of Gene and TranscriptIDs) can then be obtained for overlapping seed regions\n" +
        " \n" +
        " Note: for now, the seed regions are only matched bases on nt 2 to 8, nothing more complex is performed in\n" +
        " terms of matching.";
    }



    
    
    /**
     * This parses out the hashmap containing the run parameters for this step
     * 
     * @param configData
     * @throws Exception 
     */
    @Override
    public void parseConfigurationData(HashMap configData) throws Exception{

        logger.info(STEP_ID_STRING + ": verify configuration data");
        
        if(configData.get(ID_MINCOUNTS)==null) {
            logger.error("<" + configData.get(ID_MINCOUNTS) + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + configData.get(ID_MINCOUNTS) + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_THREADS)==null) {
            logger.error("<" + configData.get(ID_THREADS) + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + configData.get(ID_THREADS) + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_MIRBASE_VERSION)==null) {
            logger.error("<" + ID_MIRBASE_VERSION + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_MIRBASE_VERSION + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_REF_GENOME)==null) {
            logger.error("<" + ID_REF_GENOME + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_REF_GENOME + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_QUERY_GENOME)==null) {
            logger.error("<" + ID_QUERY_GENOME + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_QUERY_GENOME + "> : Missing Definition in Configuration File");
        }
                
        try{
            this.setMinCounts((Integer)configData.get(ID_MINCOUNTS));
        }
        catch(NumberFormatException exNm){
            logger.error(ID_MINCOUNTS + " <" + configData.get(ID_MINCOUNTS) + "> is not an integer");
            throw new NumberFormatException(ID_MINCOUNTS + " <" + configData.get(ID_MINCOUNTS) + "> is not an integer");
        }
        if (this.getMinCounts() <= 0){
            logger.error(ID_MINCOUNTS + " <" + configData.get(ID_MINCOUNTS) + "> must be positive");
            throw new IllegalArgumentException(ID_MINCOUNTS + " <" + configData.get(ID_MINCOUNTS) + "> must be positive");
        }
        

        try{
            this.setMiRBaseRelease((Integer) configData.get(ID_MIRBASE_VERSION));
        }
        catch(NumberFormatException exNm){
            throw new NumberFormatException(ID_MIRBASE_VERSION + " <" + configData.get(ID_MIRBASE_VERSION) + "> is not an integer");
        }        
        if (this.getMiRBaseRelease() <= 0){
            throw new IllegalArgumentException(ID_MIRBASE_VERSION + " <" + configData.get(ID_MIRBASE_VERSION) + "> must be positive integer");
        }
        
        try{
            this.setNoOfThreads((Integer) configData.get(ID_THREADS));
        }
        catch(NumberFormatException exNm){
            throw new NumberFormatException(ID_THREADS + " <" + configData.get(ID_THREADS) + "> is not an integer");
        }        
        if (this.getNoOfThreads() <= 0){
            logger.error(ID_THREADS + " <" + configData.get(ID_THREADS) + "> must be positive");
            throw new IllegalArgumentException(ID_THREADS + " <" + configData.get(ID_THREADS) + "> must be positive");
        }

        this.setQueryGenome((String) configData.get(ID_QUERY_GENOME));
        if(this.getQueryGenome().length() !=3 ){
            throw new IllegalArgumentException(ID_QUERY_GENOME + " <" + configData.get(ID_QUERY_GENOME) + "> must be a 3 letter string");            
        }
        this.setReferenceGenome((String) configData.get(ID_REF_GENOME));
        if(this.getReferenceGenome().length() !=3 ){
            throw new IllegalArgumentException(ID_REF_GENOME + " <" + configData.get(ID_REF_GENOME) + "> must be a 3 letter string");            
        }
        
 
        logger.info("passed");
    }
    
    
    
    /**
     * unzip specified file types
     * 
     * @throws IOException 
     */
    @Override
    public void execute() throws IOException, RuntimeException{

        logger.info(STEP_ID_STRING + ": execute step");      
        this.loadQueryFeaturesAndFastAFiles();
        this.loadQueryCountData();
        
        /*
        GeneralizedSuffixTree in = new GeneralizedSuffixTree();

        String word = "cacao";
        in.put(word, 0);
        */
        
        
        String gffFileMirBase = this.cleanPath(getStepInputData().getDataLocations().getMirbaseFolder() 
                + FILESEPARATOR + this.getMiRBaseRelease() + FILESEPARATOR + this.getReferenceGenome() + ".gff3");
        String faFileMirBase = gffFileMirBase.replace("gff3", "mature.fa");
        logger.info("loading miRBase Data");
        logger.info("-- GFF file is <" + gffFileMirBase + ">");
        logger.info("-- FA  file is<" + faFileMirBase + ">");
        mirBaseSet.loadMiRBaseData(this.getReferenceGenome(), gffFileMirBase, faFileMirBase);
        logger.info("read " + mirBaseSet.getNumberOfEntries() + " entries");
        logger.info("--");
        

        logger.info("reading Target Scan Data");
        String targetScanFamily = this.cleanPath(getStepInputData().getDataLocations().getTargetscanFolder()
                + FILESEPARATOR + ReferenceDataLocations.ID_TSCAN_MIRFAMILY_FILE);
        logger.info("-- miR Family file is + <" + targetScanFamily + ">");
        tScanMirFamilies.loadConservedFamilyList(targetScanFamily);
        logger.info("read " + tScanMirFamilies.getNumberOfEntries() + " entries");
        logger.info("--");
        
        
        String targetScanPrediction = this.cleanPath(getStepInputData().getDataLocations().getTargetscanFolder()
                + FILESEPARATOR + ReferenceDataLocations.ID_TSCAN_PREDICTIONS_FILE);
        logger.info("-- Predicted Target Data file is + <" + targetScanPrediction + ">");
        tscanPredictedTargets.loadPredictedTargetInfo(targetScanPrediction);
        logger.info("read " + tscanPredictedTargets.getNumberOfEntries() + " entries");
        logger.info("--\n");
        
        
        this.parseSAMFile();


        // need a list of tabbed delimited MSY mapped smallRNAs entries
        /*
            ultimately, this should be the following steps:
                1. load the MSY fasta file
                2. load the count data
                3. filter the fasta entries by counts        
                4. map the passed MSY sequences to the HSA genome
                5. parse the SAM file to extract the reads that have a mapped seed region
                6. 
            for now, lets start from step 5
        */
        
        
        logger.info("--");
        logger.info(STEP_ID_STRING + ": completed");
        
        
        
        
    }
    
    
    
    /**
     * parse the SAM file for candidate hits
     * 
     */
    private void parseSAMFile() throws IOException{
        // read in the SAM file
        // then cycle through each entry and search for conserved family match
        this.setHostVersusQuerySAMFile( this.cleanPath(this.inFolder 
                + FILESEPARATOR + getStepInputData().getProjectID() + "." 
                + this.getQueryGenome() + "_vs_" + this.getReferenceGenome() + ".sam"));
        String SAMLine = "";
        int queryRefMatches = 0;
        int queryHits = 0;
        try(BufferedReader brSM = new BufferedReader(new FileReader(new File(this.getHostVersusQuerySAMFile())))){
            while((SAMLine=brSM.readLine())!=null){
                if(SAMEntry.isHeaderLine(SAMLine)) continue;
                SAMEntry samEntry  = new SAMEntry(SAMLine);
                if(samEntry.getFlags()==0 || samEntry.getFlags()==16){
                    String qFeatureID = samEntry.getqName().split("\\|")[0].trim();
                    String mdString = samEntry.getTagValue("MD");
                    /*
                        MD:Z:15G2T4 ok
                        MD:Z:4T1A13 not ok
                        MD:Z:0C13A7 ok (first base doesnt matter)
                    */
                    // is first char a '0'?
                    Pattern mdPattern;
                    Matcher match;
                    String querySeed="";
                    if (mdString.startsWith("0")){
                        //logger.info("-- first character is '0'");
                        mdPattern = Pattern.compile("(^(0)(\\w))|((^7-9)(\\w))");
                        match = mdPattern.matcher(mdString);
                    }else{
//                        mdPattern = Pattern.compile("^[(1)(\\d)]|[(7-9)](\\w)");
                        mdPattern = Pattern.compile("(^(1)(\\d)(\\w))|(^[7-9](\\w))");
                        match = mdPattern.matcher(mdString);
                    }
                    if (match.find()){
                        // find the entry in the QueryFeature list
                        GFFEntry queryFeature = this.queryFeatures.findEntryByID(qFeatureID);
                        querySeed = (String) queryFeature.getSequence().subSequence(1, 8);
                        
                        if(Integer.parseInt(queryFeature.getAttrValue("counts")) > this.minCounts){
                            
                            try(BufferedWriter bwSH = new BufferedWriter(new FileWriter(
                                    new File(this.cleanPath(this.inFolder 
                                            + FILESEPARATOR + getStepInputData().getProjectID() + "." 
                                            + this.getQueryGenome() + "_vs_" + this.getReferenceGenome() 
                                            + "." + queryFeature.getFeatureID().replace(":", "_") + ".summary"))))){
                                logger.info(querySeed + "|" + queryFeature.getFeatureID() + "|" + queryFeature.getAttrValue("counts"));
                                logger.info("search for seed in TargetScan miR Family..");

                                bwSH.write("+" + StringUtils.repeat("-", 40) + "+\n");
                                bwSH.write("          " + queryFeature.getFeatureID() + "\n");
                                bwSH.write("+" + StringUtils.repeat("-", 40) + "+\n\n");
                                bwSH.write("          " + queryFeature.getAttrValue("counts") + "\n");
                                bwSH.write(queryFeature.getSeqID() + "\t" + queryFeature.getStart() + "\t" 
                                        + queryFeature.getStop() + "\t" + queryFeature.getStrand() + "\n");
                                bwSH.write("seed region =\t" + querySeed + "\n");
                                bwSH.write("+" + StringUtils.repeat("-", 40) + "+\n");
                                bwSH.write("          Target Hits\n");
                                bwSH.write("+" + StringUtils.repeat("-", 40) + "+\n");
                                ArrayList<String> miRFamilyHits = this.tScanMirFamilies.findSeedHits(querySeed);
                                ArrayList<String> miRTargetHits = this.tscanPredictedTargets.getGeneTargetHits(miRFamilyHits);
                                for (String hit: miRTargetHits){
                                   bwSH.write(hit + "\n"); 
                                }
                                
                                
                                bwSH.close();
                            }
                            catch(IOException exIO){
                                logger.info("error writing summary file <" + queryFeature.getFeatureID() + ".summary" + ">");
                                logger.info(exIO);
                                throw new IOException("error writing summary file <" + queryFeature.getFeatureID() + ".summary" + ">"
                                + "\n" + "see log file for details");
                            }

                        }
                        //querySeed = 
                    }
                    
                    // now try to match to a host miRNA entry
                    /*
                    if(!querySeed.isEmpty()){
                        Boolean hit = false;
                        for(MiRNAFeature miRfeature: mirBaseSet.getMiRBaseMiRNAList()){
                            if(miRfeature.getSequence()!= null && miRfeature.getSequence().substring(1, 7).equals(querySeed)){
                                GFFEntry queryHit = queryFeatures.findEntryByID(qFeatureID);
                                if(Integer.parseInt(queryHit.getAttrValue("counts")) > this.minCounts){
                                    //logger.info(querySeed);
                                    //logger.info(querySeed + ":matched " + qFeatureID + "(" + queryHit.getAttrValue("counts") + ")" + "->" + miRfeature.getName());
                                    queryRefMatches++;
                                    if(!hit) {
                                        queryHits++;
                                        hit = true;
                                    }
                                }
                            }
                        }
                    }
                    */
                }
                
            }
            logger.info("found " + queryHits + " query smallRNAs with matching seeds");
            logger.info("and a total of " +  queryRefMatches + " seed matches");
        }
        catch(IOException exIO){
            logger.error("error parsing SAM file <" + this.getHostVersusQuerySAMFile() + ">");
            logger.error("error occurred on line " + SAMLine);
            logger.error(exIO);
            throw new IOException("error parsing SAM file <" 
                    + this.getHostVersusQuerySAMFile() + ">\nsee log file for details");
        }
        
    }
    /**
     * load Query Fasta File
     * 
     * @return
     * @throws IOException 
     */
    private int loadQueryFeaturesAndFastAFiles() throws IOException, RuntimeException{
        
        this.setQueryFastaFile(this.cleanPath(this.inFolder + FILESEPARATOR + getStepInputData().getProjectID() + QUERYFA_EXTENSION));
        String faLine = "";
        //querySeedSet
        /*
            for now we just take positions 2 to 8 and store these as the seed.
            in reality it is more complex because some families appear to be shifted by +1 nt 
            e.g. the seed sequence GGAGCUC is 2 to 8 in the miR-25 family, but at 3-9 in the human counterpart
        
                >grep -A 1 msy-187803 ../analyzeStartPositions500.15.30/sweden.fasta
                >msy-187803:5 187803:5 msy 187803:5        
                UUGGGGGUUGGAGAGAAGGUCCA
        
                >grep  msy-187803 ../analyzeStartPositions500.15.30/sweden.features.tsv 
                chr3	.	smallRNA	332185	332207	.	+	.	ID=msy-187803:5;Alias=msy-187803:5;Name=msy-187803:5
        
                >grep msy-187803 ../differentialExpression.500.15.30/sweden.merged.mirna_counts.tsv
                msy-187803:5|msy-187803:5	0	0	0	0	0	0	0	0	0	0    
        
        */
        queryFeatures = new GFFSet();
        this.setQueryFeatureFile(this.cleanPath(this.inFolder + FILESEPARATOR + getStepInputData().getProjectID() + QUERYFEATURES_EXTENSION));
        logger.info("reading Query Feature file <" + this.getQueryFeatureFile() + ">");
        int lines = queryFeatures.readGFF(this.getQueryFeatureFile());
        logger.info("completed");
        logger.info("read " + lines + " entries");
        
        int faMatchCount = 0;
        int faCount = 0;
        logger.info("--");
        logger.info("reading Query Fasta file <" + this.getQueryFastaFile() + ">");
        try(BufferedReader brFA = new BufferedReader(new FileReader(new File(this.getQueryFastaFile())))){
            while((faLine=brFA.readLine())!=null){
                faCount++;
                String headerLine = faLine;
                logger.debug(headerLine);
                if (headerLine.startsWith(">")==false){
                    logger.error("error reading Query Fasta file <" + this.getQueryFastaFile() + ">");
                    logger.error("error occurred on header line " + headerLine);
                    throw new IOException("error reading Query Fasta file <" + this.getQueryFastaFile() + ">");
                }
                String seqLine = brFA.readLine();
                String thisID = faLine.substring(1).split(" ")[0].trim();
                if(queryFeatures.findEntryByID(thisID)!=null){
                    queryFeatures.findEntryByID(thisID).setSequence(seqLine.trim());
                    faMatchCount++;
                }
            }
            brFA.close();
            logger.info("read " + faCount + "lines");
            logger.info("matched " + faMatchCount + " entries in the feature file");
            logger.info("--");
        }
        catch(IOException exIO){
            logger.error("error reading Query Fasta file <" + this.getQueryFastaFile() + ">");
            logger.error("error occurred on line " + faLine);
            logger.error(exIO);
            throw new IOException("error reading Query Fasta file <" + this.getQueryFastaFile() + ">\nsee log file for details");
        }
        return queryFeatures.getNoOfEntries();
    }
    
    
    
    
    /**
     * load Query Count data file
     * 
     * @throws IOException 
     */
    private void loadQueryCountData() throws IOException{
        int cCount = 0;
        int ctMatchCount = 0;
        int dualMatchCount = 0;
        String countLine = "";
        this.setQueryCountDataFile(this.cleanPath(this.inFolder + FILESEPARATOR + getStepInputData().getProjectID() + QUERYCOUNTS_EXTENSION));
        logger.info("reading Query count file <" + this.getQueryCountDataFile() + ">");
        try(BufferedReader brCF = new BufferedReader(new FileReader(new File(this.getQueryCountDataFile())))){
            brCF.readLine();
            while((countLine=brCF.readLine())!=null){
                cCount++;
                String counts[] = countLine.split("\t");
                int totalCounts = 0;
                for(int i=1; i< counts.length; i++){
                    totalCounts += Integer.parseInt(counts[i]);
                }
                int avgCounts = totalCounts / (counts.length-1);
                String thisFeatureID = counts[0].split("\\|")[0].trim();
                if(queryFeatures.findEntryByID(thisFeatureID)!=null){
                    queryFeatures.findEntryByID(thisFeatureID).addAttr("avgcounts", Integer.toString(avgCounts));
                    ctMatchCount++;
                    if(queryFeatures.findEntryByID(thisFeatureID).getAttr().contains("seq=")){
                        dualMatchCount++;
                    }
                }
            }
            brCF.close();
            logger.info("read " + cCount + " lines");
            logger.info("matched " + ctMatchCount + " entries in the feature file");
            logger.info("matched " + dualMatchCount + " entries in the feature file with 'seq=' attr");
            logger.info("--");
            
        }
        catch(IOException exIO){
            logger.error("error reading Query count file <" + this.getQueryCountDataFile() + ">");
            logger.error("error occurred on line " + countLine);
            logger.error(exIO);
            throw new IOException("error reading Query Fasta file <" + this.getQueryCountDataFile() + ">\nsee log file for details");
        }

        
    }
    /**
     * this should be called prior to executing the step.
     * check unzip software exists and input files are available
     * 
     * @throws IOException
     */
    @Override
    public void verifyInputData() throws IOException{
        
        logger.info("verify input data");  

        this.setPaths();
        this.setQueryFeatureFile(this.cleanPath(this.inFolder + FILESEPARATOR + getStepInputData().getProjectID() + QUERYFEATURES_EXTENSION));        
        if(new File(this.getQueryFeatureFile()).exists() == false){
            throw new IOException(STEP_ID_STRING + ": Query feature file < " + this.getQueryFeatureFile() +"> not found");
        }

        this.setQueryFastaFile(this.cleanPath(this.inFolder + FILESEPARATOR + getStepInputData().getProjectID() + QUERYFA_EXTENSION));        
        if(new File(this.getQueryFastaFile()).exists() == false){
            throw new IOException(STEP_ID_STRING + ": Query fasta file < " + this.getQueryFastaFile() +"> not found");
        }
        
        this.setQueryCountDataFile(this.cleanPath(this.inFolder + FILESEPARATOR + getStepInputData().getProjectID() + QUERYCOUNTS_EXTENSION));
        if(new File(this.getQueryCountDataFile()).exists() == false){
            throw new IOException(STEP_ID_STRING + ": Query counts file < " + this.getQueryCountDataFile() +"> not found");
        }
        
        
        this.setHostVersusQuerySAMFile( this.cleanPath(this.inFolder 
                + FILESEPARATOR + getStepInputData().getProjectID() + "." 
                + this.getQueryGenome() + "_vs_" + this.getReferenceGenome() + ".sam"));
        if(new File(this.getHostVersusQuerySAMFile()).exists() == false){
            throw new IOException(STEP_ID_STRING + ": SAM hits file < " + this.getHostVersusQuerySAMFile() +"> not found");
        }
        
        /*
        target scan files
        */    
                
        String gffFileMirBase = this.cleanPath(getStepInputData().getDataLocations().getMirbaseFolder() 
                + FILESEPARATOR + this.getMiRBaseRelease() + FILESEPARATOR + this.getReferenceGenome() + ".gff3");        
        if(new File(gffFileMirBase).exists() == false){
            throw new IOException(STEP_ID_STRING + ": miRBase GFF file not found at location < " + gffFileMirBase +">");
        }
        
        String faFileMirBase = gffFileMirBase.replace("gff3", "mature.fa");
        if(new File(faFileMirBase).exists() == false){
            throw new IOException(STEP_ID_STRING + ": miRBase Fasta file not found at location < " + faFileMirBase +">");
        }
    
        String targetScanFamily = this.cleanPath(getStepInputData().getDataLocations().getTargetscanFolder()
                + FILESEPARATOR + ReferenceDataLocations.ID_TSCAN_MIRFAMILY_FILE);
        if(new File(targetScanFamily).exists() == false){
            throw new IOException(STEP_ID_STRING + ": TargetScan miR Family file not found at location < " + faFileMirBase +">");
        }
        
        
        String targetScanPrediction = this.cleanPath(getStepInputData().getDataLocations().getTargetscanFolder()
                + FILESEPARATOR + ReferenceDataLocations.ID_TSCAN_MIRFAMILY_FILE);
        if(new File(targetScanPrediction).exists() == false){
            throw new IOException(STEP_ID_STRING + ": TargetScan Predicted Targets file not found at location < " + faFileMirBase +">");
        }
            
                            
        // check the data files are present
                        
            


    }
    
    
    
    @Override
    public void verifyOutputData(){
        logger.info("no output verification required");
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
        
        HashMap<String, Object> configData = new HashMap();
        
        configData.put(ID_REF_GENOME, "hsa");
        configData.put(ID_QUERY_GENOME, "hsa");
        configData.put(ID_MIRBASE_VERSION, 20);
        configData.put(ID_MINCOUNTS, 1000);
        configData.put(ID_THREADS, 4);
        
        return configData;
        
    }

    
    
    
    
    
    /**
     * @return the noOfThreads
     */
    public int getNoOfThreads() {
        return noOfThreads;
    }

    /**
     * @param noOfThreads the noOfThreads to set
     */
    public void setNoOfThreads(int noOfThreads) {
        this.noOfThreads = noOfThreads;
    }

    /**
     * @return the unzipSoftware
     */
    public String getUnzipSoftware() {
        return unzipSoftware;
    }

    /**
     * @param unzipSoftware the unzipSoftware to set
     */
    public void setUnzipSoftware(String unzipSoftware) {
        this.unzipSoftware = unzipSoftware;
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
     * @return the referenceGenome
     */
    public String getReferenceGenome() {
        return referenceGenome;
    }

    /**
     * @param referenceGenome the referenceGenome to set
     */
    public void setReferenceGenome(String referenceGenome) {
        this.referenceGenome = referenceGenome;
    }

    /**
     * @return the queryGenome
     */
    public String getQueryGenome() {
        return queryGenome;
    }

    /**
     * @param queryGenome the queryGenome to set
     */
    public void setQueryGenome(String queryGenome) {
        this.queryGenome = queryGenome;
    }

    /**
     * @return the queryCountDataFile
     */
    public String getQueryCountDataFile() {
        return queryCountDataFile;
    }

    /**
     * @param queryCountDataFile the queryCountDataFile to set
     */
    public void setQueryCountDataFile(String queryCountDataFile) {
        this.queryCountDataFile = queryCountDataFile;
    }

    /**
     * @return the queryFastaFile
     */
    public String getQueryFastaFile() {
        return queryFastaFile;
    }

    /**
     * @param queryFastaFile the queryFastaFile to set
     */
    public void setQueryFastaFile(String queryFastaFile) {
        this.queryFastaFile = queryFastaFile;
    }

    /**
     * @return the hostVersusQuerySAMFile
     */
    public String getHostVersusQuerySAMFile() {
        return hostVersusQuerySAMFile;
    }

    /**
     * @param hostVersusQuerySAMFile the hostVersusQuerySAMFile to set
     */
    public void setHostVersusQuerySAMFile(String hostVersusQuerySAMFile) {
        this.hostVersusQuerySAMFile = hostVersusQuerySAMFile;
    }

    /**
     * @return the queryFeatureFile
     */
    public String getQueryFeatureFile() {
        return queryFeatureFile;
    }

    /**
     * @param queryFeatureFile the queryFeatureFile to set
     */
    public void setQueryFeatureFile(String queryFeatureFile) {
        this.queryFeatureFile = queryFeatureFile;
    }

    /**
     * @return the minCounts
     */
    public int getMinCounts() {
        return minCounts;
    }

    /**
     * @param minCounts the minCounts to set
     */
    public void setMinCounts(int minCounts) {
        this.minCounts = minCounts;
    }

   
}
